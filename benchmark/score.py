#!/usr/bin/env python3
"""Score the benchmark: does the pipeline send each planted problem to review and let clean runs pass?

    python benchmark/score.py                                   # gates that need no model
    python benchmark/score.py --ask --region us-east-1          # also ask Claude on Bedrock, then score
    python benchmark/score.py --answers benchmark/answers/<dir> # score answers saved by an earlier --ask

Without a model, only the trust and fact gates can be scored; scenarios that need a judgment are
listed as pending. With answers, every scenario gets a final verdict.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / "benchmark"
QUESTIONS = ROOT / "skills" / "llomics-limma" / "references" / "questions.contrast.json"


def load_judge():
    path = ROOT / "skills" / "llomics-judge" / "scripts" / "judge_bedrock.py"
    spec = importlib.util.spec_from_file_location("judge_bedrock", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gate_only(jb, payload: dict[str, Any], gates: dict[str, Any]) -> list[str]:
    unit = payload["units"][0]
    return jb.trust_reasons(payload, unit, gates) + jb.fact_reasons(unit, gates)


def layer(reason: str) -> str:
    return "trust" if reason.startswith("trust:") else "fact"


def verdict_without_model(scenario: dict[str, Any], reasons: list[str]) -> str:
    expected, caught_by = scenario["expected"], scenario["caught_by"]
    if expected == "auto":
        return "FALSE ALARM" if reasons else "ok so far (model decides)"
    if reasons:
        layers = sorted({layer(r) for r in reasons})
        return "caught" if caught_by in layers else f"caught by {'+'.join(layers)} (expected {caught_by})"
    if caught_by == "model":
        return "pending (needs the model)"
    if caught_by == "none":
        return "not caught (known blind spot)"
    return "MISSED"


def verdict_with_model(scenario: dict[str, Any], status: str) -> str:
    if status == scenario["expected"]:
        return "correct"
    if scenario["expected"] == "auto":
        return "FALSE ALARM"
    return "MISSED (known blind spot)" if scenario["caught_by"] == "none" else "MISSED"


def evaluate(gates: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    """Gate-only results for every scenario. Used by the tests as well as the command line."""
    jb = load_judge()
    cfg = json.loads(QUESTIONS.read_text(encoding="utf-8"))
    gates = gates or cfg["gates"]
    scenarios = json.loads((BENCH / "scenarios.json").read_text(encoding="utf-8"))["scenarios"]
    rows = []
    for s in scenarios:
        payload = json.loads((BENCH / "payloads" / f"{s['id']}.llomics.json").read_text(encoding="utf-8"))
        reasons = gate_only(jb, payload, gates)
        rows.append({**s, "reasons": reasons, "without_model": verdict_without_model(s, reasons)})
    return rows


def ask(args: argparse.Namespace, out_dir: Path) -> None:
    jb = load_judge()
    cfg = json.loads(QUESTIONS.read_text(encoding="utf-8"))
    call = jb.make_caller(args)
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in sorted((BENCH / "payloads").glob("*.llomics.json")):
        raw = path.read_bytes()
        payload = json.loads(raw)
        results = jb.judge_units(payload, cfg, call, args.runs, args.workers)
        record = {"payload_file": path.name, "payload_sha256": hashlib.sha256(raw).hexdigest(),
                  "model_requested": args.model, "client": args.client, "runs_per_unit": args.runs,
                  "gates": cfg["gates"], "units": results}
        target = out_dir / path.name.replace(".llomics.json", ".answers.json")
        target.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"asked {path.name}: {results[0]['status']} {results[0]['reasons']}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--answers", type=Path, help="directory of <id>.answers.json files to score")
    parser.add_argument("--ask", action="store_true", help="call Claude on Bedrock for every scenario first")
    parser.add_argument("--region", default=os.environ.get("AWS_REGION") or os.environ.get("AWS_DEFAULT_REGION"))
    parser.add_argument("--model", default=os.environ.get("LLOMICS_MODEL", "anthropic.claude-opus-5-5"))
    parser.add_argument("--client", choices=["mantle", "legacy"], default="mantle")
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max-tokens", type=int, default=16000)
    args = parser.parse_args()

    answers_dir = args.answers
    if args.ask:
        if not args.region:
            print("ERROR  set --region or AWS_REGION")
            return 1
        answers_dir = answers_dir or BENCH / "answers" / re.sub(r"[^A-Za-z0-9._-]", "_", args.model)
        try:
            ask(args, answers_dir)
        except ImportError:
            print("ERROR  the anthropic SDK is not installed: pip install anthropic")
            return 1

    rows = evaluate()
    print(f"\n{'scenario':<22} {'expected':<8} {'should be caught by':<20} without a model")
    print("-" * 100)
    for r in rows:
        detail = f"  [{', '.join(r['reasons'])}]" if r["reasons"] else ""
        print(f"{r['id']:<22} {r['expected']:<8} {str(r['caught_by'] or '-'):<20} {r['without_model']}{detail}")

    bad = [r for r in rows if r["without_model"] in ("MISSED", "FALSE ALARM")]
    pending = [r for r in rows if r["without_model"].startswith("pending")]
    print(f"\nGates alone: {len(bad)} wrong, {len(pending)} pending on the model, "
          f"{sum(r['without_model'].startswith('not caught') for r in rows)} known blind spot(s).")

    if answers_dir:
        print(f"\n{'scenario':<22} {'expected':<8} {'final':<8} {'verdict':<26} dominant_signal (share) / reasons")
        print("-" * 110)
        n_correct = 0
        for r in rows:
            path = answers_dir / f"{r['id']}.answers.json"
            if not path.exists():
                print(f"{r['id']:<22} {r['expected']:<8} {'-':<8} no answers file")
                continue
            unit = json.loads(path.read_text(encoding="utf-8"))["units"][0]
            verdict = verdict_with_model(r, unit["status"])
            n_correct += verdict == "correct"
            signal = unit["summary"].get("dominant_signal", {})
            shown = f"{signal.get('top')} ({signal.get('share')})" if signal else "-"
            print(f"{r['id']:<22} {r['expected']:<8} {unit['status']:<8} {verdict:<26} {shown} {unit['reasons']}")
        print(f"\nWith the model: {n_correct} of {len(rows)} scenarios decided correctly.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
