#!/usr/bin/env python3
"""Ask Claude on Amazon Bedrock typed questions about each unit of an llomics payload, then gate in code.

Usage:
    python judge_bedrock.py payload.json questions.json --out answers.json --region us-east-1
    python judge_bedrock.py payload.json questions.json --dry-run          # no model call, no SDK needed

Requires `pip install anthropic` for real runs. Credentials come from the standard AWS chain
(environment variables, shared profile, or an attached role); none are read or written here.

Each unit is asked `--runs` times. The share of runs that agree stands in for a probability,
because a general LLM does not return calibrated probabilities. Gates are evaluated in this order:
trust checks, then facts from the payload, then model judgments. Any failed gate sends the unit to
review with a reason code.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import operator
import os
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable

DEFAULT_MODEL = "anthropic.claude-opus-5-5"

DEFAULT_SYSTEM = (
    "You are helping a statistician review the output of an omics analysis. The user message is "
    "JSON with three parts. `context` holds the scientific question, the filters that were applied, "
    "the dataset size, and the status of the quality checks. `unit` is one unit of the analysis, "
    "for example one contrast with its hit counts, its top features, and a panel of features the "
    "analyst expected to matter. `questions` lists what to answer.\n\n"
    "Answer each question about this unit only, from the state provided and your knowledge of "
    "biology. In feature rows, `effect` is the effect size (log2 fold change for a differential "
    "analysis), `stat` is the test statistic, and `p_adj` is the adjusted p-value. In `panel`, a "
    "null value means the feature was not measured or was filtered out; it does not mean zero. "
    "The top-feature lists are a ranked sample: `tail` says how many features were eligible.\n\n"
    "Hit counts, thresholds, and quality checks are evaluated separately in code. Your answers "
    "should reflect what the features are, not how many there are or how significant they look."
)

OPS: dict[str, Callable[[Any, Any], bool]] = {
    "<": operator.lt, "<=": operator.le, ">": operator.gt, ">=": operator.ge, "==": operator.eq,
}
VOTED_SHAPES = {"bool", "choice", "score"}


def build_schema(questions: list[dict[str, Any]]) -> dict[str, Any]:
    """JSON schema with one required property per question."""
    props: dict[str, Any] = {}
    for q in questions:
        shape = q["shape"]
        if shape == "bool":
            props[q["id"]] = {"type": "boolean", "description": q["text"]}
        elif shape == "choice":
            props[q["id"]] = {"type": "string", "enum": list(q["options"]), "description": q["text"]}
        elif shape == "score":
            props[q["id"]] = {"type": "string", "enum": list(q["levels"]), "description": q["text"]}
        elif shape == "text":
            props[q["id"]] = {"type": "string", "description": q["text"]}
        else:
            raise ValueError(f"question {q['id']}: unknown shape '{shape}'")
    return {"type": "object", "properties": props, "required": list(props), "additionalProperties": False}


def build_request(payload: dict[str, Any], unit: dict[str, Any], cfg: dict[str, Any]) -> dict[str, Any]:
    """The exact request for one unit: system prompt, user content, and output schema."""
    run = payload["run"]
    state = {
        "context": {
            "method": payload["method"],
            "question": run.get("question"),
            "tool": run.get("tool"),
            "steps": run.get("steps"),
            "filters": run.get("filters"),
            "dataset": payload["dataset"],
            "quality_checks": payload["trust"],
        },
        "unit": unit,
        "questions": [
            {k: q[k] for k in ("id", "text", "options", "levels") if k in q} for q in cfg["questions"]
        ],
    }
    return {
        "system": cfg.get("system", DEFAULT_SYSTEM),
        "user": json.dumps(state, indent=1, ensure_ascii=False),
        "schema": build_schema(cfg["questions"]),
    }


def parse_response(response: Any) -> dict[str, Any]:
    """Turn a Messages API response into {ok, answers} or {ok: False, error}. A failure is never an answer."""
    stop = response.stop_reason
    if stop == "refusal":
        details = getattr(response, "stop_details", None)
        category = getattr(details, "category", None) if details else None
        return {"ok": False, "error": f"refusal ({category})" if category else "refusal"}
    if stop == "max_tokens":
        return {"ok": False, "error": "truncated at max_tokens"}
    text = next((b.text for b in response.content if b.type == "text"), None)
    if text is None:
        return {"ok": False, "error": f"no text block (stop_reason={stop})"}
    try:
        answers = json.loads(text)
    except json.JSONDecodeError as exc:
        return {"ok": False, "error": f"invalid JSON: {exc}"}
    return {"ok": True, "answers": answers, "model": response.model}


def make_caller(args: argparse.Namespace) -> Callable[[dict[str, Any]], dict[str, Any]]:
    """Build the function that sends one request to Bedrock. Imports the SDK only when needed."""
    import anthropic

    if args.client == "mantle":
        from anthropic import AnthropicBedrockMantle

        client = AnthropicBedrockMantle(aws_region=args.region)
    else:  # legacy bedrock-runtime InvokeModel path
        from anthropic import AnthropicBedrock

        client = AnthropicBedrock(aws_region=args.region)

    def call(request: dict[str, Any]) -> dict[str, Any]:
        try:
            response = client.messages.create(
                model=args.model,
                max_tokens=args.max_tokens,
                system=request["system"],
                messages=[{"role": "user", "content": request["user"]}],
                output_config={"format": {"type": "json_schema", "schema": request["schema"]}},
            )
        except anthropic.RateLimitError as exc:
            return {"ok": False, "error": f"rate_limited: {exc}"}
        except anthropic.APIStatusError as exc:
            return {"ok": False, "error": f"http_{exc.status_code}: {exc.message}"}
        except anthropic.APIConnectionError as exc:
            return {"ok": False, "error": f"connection: {exc}"}
        return parse_response(response)

    return call


def summarise(runs: list[dict[str, Any]], questions: list[dict[str, Any]]) -> dict[str, Any]:
    """Per question: the answers given, the most common one, and the share of successful runs that gave it."""
    ok = [r["answers"] for r in runs if r["ok"]]
    out: dict[str, Any] = {}
    for q in questions:
        values = [a.get(q["id"]) for a in ok]
        if q["shape"] not in VOTED_SHAPES:
            out[q["id"]] = {"answers": values}
            continue
        counts = Counter(json.dumps(v) for v in values)
        ranked = counts.most_common()
        top = json.loads(ranked[0][0]) if ranked else None
        share = ranked[0][1] / len(values) if values else 0.0
        second = ranked[1][1] / len(values) if len(ranked) > 1 else 0.0
        out[q["id"]] = {
            "top": top,
            "share": round(share, 3),
            "margin": round(share - second, 3),
            # keys are the answers: option text as is, booleans as "true"/"false"
            "counts": {(json.loads(k) if k.startswith('"') else k): n for k, n in ranked},
        }
    return out


def trust_reasons(payload: dict[str, Any], unit: dict[str, Any], gates: dict[str, Any]) -> list[str]:
    """Reason codes for blocking trust checks.

    ``blocking_checks`` is a list of check names (each blocks on ``blocking_statuses``), or an object
    mapping each check name to the statuses that block it, so that a warning can block for one check
    and not for another.
    """
    default = gates.get("blocking_statuses", ["fail", "not_run"])
    spec = gates.get("blocking_checks", [])
    per_check = spec if isinstance(spec, dict) else {name: default for name in spec}
    reasons = []
    for name, blocking in per_check.items():
        statuses = [c["status"] for c in payload["trust"] + unit.get("trust", []) if c["check"] == name]
        if not statuses:
            reasons.append(f"trust:{name}:absent")
        for status in sorted(set(statuses) & set(blocking)):
            reasons.append(f"trust:{name}:{status}")
    return reasons


def fact_reasons(unit: dict[str, Any], gates: dict[str, Any]) -> list[str]:
    reasons = []
    for gate in gates.get("facts", []):
        value = unit.get("facts", {}).get(gate["fact"])
        if value is None:
            reasons.append(f"fact:{gate['fact']}:missing")
        elif OPS[gate["op"]](value, gate["value"]):
            reasons.append(gate["reason"])
    return reasons


def judgment_reasons(summary: dict[str, Any], runs: list[dict[str, Any]], gates: dict[str, Any]) -> list[str]:
    reasons = []
    n_ok = sum(r["ok"] for r in runs)
    if n_ok < len(runs):
        reasons.append(f"no_answer:{len(runs) - n_ok}_of_{len(runs)}_runs")
    if n_ok == 0:
        return reasons
    for gate in gates.get("judgments", []):
        s = summary[gate["question"]]
        accepted = s["top"] in gate["accept"]
        if not accepted or s["share"] < gate.get("min_share", 0.8):
            reasons.append(gate["reason"])
    return reasons


def judge_units(
    payload: dict[str, Any],
    cfg: dict[str, Any],
    call: Callable[[dict[str, Any]], dict[str, Any]] | None,
    n_runs: int,
    workers: int,
    limit: int | None = None,
) -> list[dict[str, Any]]:
    """Ask every unit n_runs times (unless call is None) and apply the gates."""
    gates = cfg.get("gates", {})
    units = payload["units"][:limit] if limit else payload["units"]
    requests = [build_request(payload, u, cfg) for u in units]

    all_runs: list[list[dict[str, Any]]] = [[] for _ in units]
    if call is not None:
        jobs = [(i, requests[i]) for i in range(len(units)) for _ in range(n_runs)]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for (i, _), result in zip(jobs, pool.map(lambda job: call(job[1]), jobs)):
                all_runs[i].append(result)

    results = []
    for unit, request, runs in zip(units, requests, all_runs):
        reasons = trust_reasons(payload, unit, gates) + fact_reasons(unit, gates)
        summary = summarise(runs, cfg["questions"]) if runs else {}
        if call is not None:
            reasons += judgment_reasons(summary, runs, gates)
        status = "not_asked" if call is None and not reasons else ("review" if reasons else "auto")
        results.append({
            "unit_id": unit["unit_id"],
            "status": status,
            "reasons": reasons,
            "summary": summary,
            "runs": runs,
            "request": request,
        })
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("payload", type=Path)
    parser.add_argument("questions", type=Path)
    parser.add_argument("--out", type=Path, default=Path("answers.json"))
    parser.add_argument("--region", default=os.environ.get("AWS_REGION") or os.environ.get("AWS_DEFAULT_REGION"))
    parser.add_argument("--model", default=os.environ.get("LLOMICS_MODEL", DEFAULT_MODEL),
                        help=f"Bedrock model ID enabled in your account (default: {DEFAULT_MODEL})")
    parser.add_argument("--client", choices=["mantle", "legacy"], default="mantle",
                        help="mantle: Bedrock Messages-API endpoint (default). legacy: bedrock-runtime InvokeModel")
    parser.add_argument("--runs", type=int, default=3, help="requests per unit; the vote share needs more than one")
    parser.add_argument("--workers", type=int, default=4, help="concurrent requests")
    parser.add_argument("--max-tokens", type=int, default=16000)
    parser.add_argument("--limit", type=int, help="only the first N units (smoke test)")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the first request and the trust/fact gates; no model call")
    args = parser.parse_args()

    raw = args.payload.read_bytes()
    payload = json.loads(raw)
    cfg = json.loads(args.questions.read_text(encoding="utf-8"))
    expected_type = cfg.get("unit_type")

    if expected_type:
        mismatched = [u["unit_id"] for u in payload["units"] if u.get("unit_type") != expected_type]
        if mismatched:
            print(f"ERROR  questions are for unit_type '{expected_type}'; other units present: {mismatched}")
            return 1

    if args.dry_run:
        results = judge_units(payload, cfg, None, 0, 1, args.limit)
        first = results[0]["request"]
        print("=== system ===\n" + first["system"])
        print("\n=== user (unit " + str(results[0]["unit_id"]) + ") ===\n" + first["user"])
        print("\n=== output schema ===\n" + json.dumps(first["schema"], indent=1))
        print("\n=== gates that do not need the model ===")
        for r in results:
            print(f"{r['unit_id']}: {', '.join(r['reasons']) if r['reasons'] else 'no trust or fact gate failed'}")
        chars = [len(r["request"]["system"]) + len(r["request"]["user"]) for r in results]
        print(f"\n{len(results)} unit(s); request size {min(chars)}-{max(chars)} characters. Nothing was sent.")
        return 0

    if not args.region:
        print("ERROR  set --region or AWS_REGION")
        return 1
    try:
        call = make_caller(args)
    except ImportError:
        print("ERROR  the anthropic SDK is not installed: pip install anthropic")
        return 1

    results = judge_units(payload, cfg, call, args.runs, args.workers, args.limit)

    record = {
        "llomics_answers": "0.1",
        "payload_file": str(args.payload),
        "payload_sha256": hashlib.sha256(raw).hexdigest(),
        "questions_file": str(args.questions),
        "model_requested": args.model,
        "client": args.client,
        "region": args.region,
        "runs_per_unit": args.runs,
        "gates": cfg.get("gates", {}),
        "units": results,
    }
    args.out.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    for r in results:
        flag = "AUTO  " if r["status"] == "auto" else "REVIEW"
        print(f"{flag} {r['unit_id']}: {', '.join(r['reasons']) if r['reasons'] else 'all gates passed'}")
    n_errors = sum(1 for r in results for run in r["runs"] if not run["ok"])
    print(f"Wrote {args.out}: {len(results)} unit(s), {n_errors} failed request(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
