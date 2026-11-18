"""Glue between the web demo and llomics: JSON in, JSON out.

The page runs this file in the browser with Pyodide, next to unmodified copies of the Python
package, the judge script and the benchmark, laid out as in the repository. It runs under CPython
too, which is how tests/test_web.py checks it. Nothing here decides what goes into a payload: that
stays in llomics.core. This file reads tables and form fields, calls the package, and shapes the
results for display.
"""

from __future__ import annotations

import csv
import importlib.util
import io
import json
import math
import random
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT / "python") not in sys.path:
    sys.path.insert(0, str(ROOT / "python"))

import llomics  # noqa: E402
from llomics import core, enrich  # noqa: E402
from llomics import validate as validate_payload  # noqa: E402

QUESTIONS = ROOT / "skills" / "llomics-limma" / "references" / "questions.contrast.json"
ANSWERS = ROOT / "benchmark" / "answers"
MAX_TABLE_ROWS = 200_000


def _module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


JUDGE = _module("judge_bedrock", ROOT / "skills" / "llomics-judge" / "scripts" / "judge_bedrock.py")
SCORE = _module("llomics_benchmark_score", ROOT / "benchmark" / "score.py")


def _questions() -> dict[str, Any]:
    return json.loads(QUESTIONS.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------------------------
# Entry point for the page
# ---------------------------------------------------------------------------------------------

def call(name: str, args_json: str = "{}") -> str:
    """Run one of the public functions below. Returns {"ok": true, "result": ...} or {"ok": false, "error": ...}."""
    fn = PUBLIC.get(name)
    if fn is None:
        return json.dumps({"ok": False, "error": f"unknown function: {name}"})
    try:
        result = fn(**json.loads(args_json or "{}"))
        return json.dumps({"ok": True, "result": result}, allow_nan=False, default=str)
    except (ValueError, KeyError, TypeError) as exc:
        message = exc.args[0] if isinstance(exc, KeyError) and exc.args else str(exc)
        return json.dumps({"ok": False, "error": str(message)})


def info() -> dict[str, Any]:
    return {
        "llomics": llomics.__version__,
        "contract": llomics.VERSION,
        "python": ".".join(str(x) for x in sys.version_info[:3]),
        "answers": _answer_sets(),
    }


# ---------------------------------------------------------------------------------------------
# Planted problems: the benchmark, scored the same way as benchmark/score.py
# ---------------------------------------------------------------------------------------------

def _answer_sets() -> list[str]:
    if not ANSWERS.is_dir():
        return []
    return sorted((p.name for p in ANSWERS.iterdir() if p.is_dir() and any(p.glob("*.answers.json"))), reverse=True)


def _labels(rows: list[dict[str, Any]], n: int = 8) -> list[str]:
    return [r.get("symbol") or r["feature"] for r in rows[:n]]


SEVERITY = ["fail", "warn", "not_run", "pass", "not_applicable"]


def _worst_per_check(checks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One entry per check name, in first-seen order, with the most severe status among run and unit level."""
    out: dict[str, dict[str, Any]] = {}
    for c in checks:
        rank = SEVERITY.index(c["status"]) if c["status"] in SEVERITY else 0
        seen = out.get(c["check"])
        if seen is None or rank < SEVERITY.index(seen["status"]):
            out[c["check"]] = {"check": c["check"], "status": c["status"], "note": c.get("note")}
    return list(out.values())


def gallery(answers: str | None = None) -> dict[str, Any]:
    """Every benchmark scenario with what the gates decide, and the model's saved answers if any.

    Saved runs are gated again here with the current gate settings, so the page shows what the code
    in this version decides, not what an older version decided.
    """
    cfg = _questions()
    sets = _answer_sets()
    answer_set = answers if answers in sets else (sets[0] if sets else None)
    rows = SCORE.evaluate()
    scenarios = []
    for row in rows:
        path = ROOT / "benchmark" / "payloads" / f"{row['id']}.llomics.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        unit = payload["units"][0]
        facts = unit.get("facts", {})
        request = JUDGE.build_request(payload, unit, cfg)
        item = {
            "id": row["id"],
            "truth": row["truth"],
            "expected": row["expected"],
            "caught_by": row["caught_by"],
            "reasons": row["reasons"],
            "without_model": row["without_model"],
            "question": payload["run"].get("question"),
            "definition": unit.get("definition"),
            "n_hits": facts.get("n_hits"),
            "n_up": facts.get("n_up"),
            "n_down": facts.get("n_down"),
            "top_up": _labels(unit["top_features"]["up"]),
            "top_down": _labels(unit["top_features"]["down"]),
            "technical_flags": facts.get("technical_flags", []),
            "trust": _worst_per_check(payload["trust"] + unit.get("trust", [])),
            "payload_path": f"benchmark/payloads/{row['id']}.llomics.json",
            "payload_chars": len(json.dumps(payload)),
            "request_chars": len(request["system"]) + len(request["user"]),
            "model": None,
        }
        if answer_set:
            item["model"] = _replay(payload, cfg, row, ROOT / "benchmark" / "answers" / answer_set)
        scenarios.append(item)

    problems = [s for s in scenarios if s["expected"] == "review"]
    clean = [s for s in scenarios if s["expected"] == "auto"]
    summary = {
        "n_scenarios": len(scenarios),
        "n_problems": len(problems),
        "n_clean": len(clean),
        "caught_by_gates": sum(1 for s in problems if s["reasons"]),
        "false_alarms": sum(1 for s in clean if s["reasons"]),
        "pending_model": sum(1 for s in scenarios if s["without_model"].startswith("pending")),
        "blind_spots": sum(1 for s in scenarios if s["without_model"].startswith("not caught")),
        "answers": answer_set,
    }
    replayed = [s for s in scenarios if s["model"]]
    if replayed:
        summary["with_model_correct"] = sum(1 for s in replayed if s["model"]["verdict"] == "correct")
        summary["with_model_total"] = len(replayed)
    return {"summary": summary, "scenarios": scenarios}


def _replay(payload: dict[str, Any], cfg: dict[str, Any], row: dict[str, Any], folder: Path) -> dict[str, Any] | None:
    path = folder / f"{row['id']}.answers.json"
    if not path.exists():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    runs = [u.get("runs", []) for u in record["units"]]
    result = JUDGE.gate_units(payload, cfg, runs, limit=len(runs))[0]
    signal = result["summary"].get("dominant_signal") or {}
    models = sorted({r.get("model") for r in result["runs"] if r.get("ok") and r.get("model")})
    return {
        "set": folder.name,
        "model": ", ".join(models) or record.get("model_requested"),
        "n_runs": len(result["runs"]),
        "status": result["status"],
        "reasons": result["reasons"],
        "verdict": SCORE.verdict_with_model(row, result["status"]),
        "dominant_signal": {"top": signal.get("top"), "share": signal.get("share")},
        "summary": result["summary"],
    }


# ---------------------------------------------------------------------------------------------
# Inspect any payload: validator, units, gates, and the request that would be sent
# ---------------------------------------------------------------------------------------------

def inspect(text: str) -> dict[str, Any]:
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        return {"parse_error": f"not valid JSON: {exc}"}
    errors, warnings = validate_payload(payload)
    out: dict[str, Any] = {"errors": errors, "warnings": warnings, "chars": len(json.dumps(payload))}
    if errors or not isinstance(payload, dict):
        return out

    cfg = _questions()
    unit_type = cfg.get("unit_type")
    units = payload["units"]
    askable = bool(units) and all(u.get("unit_type") == unit_type for u in units)
    gated = {r["unit_id"]: r for r in JUDGE.gate_units(payload, cfg, None)} if askable else {}

    out.update({
        "llomics": payload["llomics"],
        "method": payload["method"],
        "run": payload["run"],
        "dataset": payload["dataset"],
        "trust": payload["trust"],
        "askable": askable,
        "questions_note": None if askable else (
            f"The shipped questions file is for {unit_type} units, and this payload has "
            f"{' and '.join(sorted({str(u.get('unit_type')) for u in units}))} units, so the page shows the units "
            "without a request or gates. A questions file for cluster units is not written yet."),
        "units": [],
    })
    for unit in units:
        g = gated.get(unit["unit_id"])
        entry = {"unit": unit, "chars": len(json.dumps(unit)), "gate": None, "request": None}
        if g:
            entry["gate"] = {"status": g["status"], "reasons": g["reasons"]}
            entry["request"] = g["request"]
        out["units"].append(entry)
    return out


# ---------------------------------------------------------------------------------------------
# Build a payload from a differential table
# ---------------------------------------------------------------------------------------------

# Column names used by common tools, mapped onto the standard columns (compared in lower case).
ALIASES: dict[str, tuple[str, ...]] = {
    "feature": ("feature", "gene", "genes", "gene_id", "geneid", "gene_name", "id", "names", "ensembl",
                "ensembl_id", "ensembl_gene_id", "row", "rownames", "unnamed: 0", ""),
    "symbol": ("symbol", "gene_symbol", "hgnc_symbol", "external_gene_name", "mgi_symbol"),
    "effect": ("effect", "logfc", "log2foldchange", "log2fc", "lfc", "avg_log2fc", "avg_logfc", "logfoldchanges"),
    "p": ("p", "pvalue", "p.value", "pval", "p_val", "pvals", "p_value", "p-value"),
    "p_adj": ("p_adj", "padj", "adj.p.val", "fdr", "p_val_adj", "pvals_adj", "qvalue", "q_value", "q.value",
              "adj_pval", "adj.pval", "p.adjust", "p_adjusted", "bh"),
    "stat": ("stat", "t", "scores", "z", "zscore", "z_score", "wald", "lr"),
    "ave_expr": ("ave_expr", "aveexpr", "logcpm", "mean_expr", "avg_expr"),
}
TOOL_HINTS = {
    "limma": {"logfc", "adj.p.val", "aveexpr"},
    "DESeq2": {"log2foldchange", "padj", "basemean"},
    "edgeR": {"logfc", "logcpm", "fdr"},
    "Seurat": {"avg_log2fc", "p_val_adj"},
    "scanpy": {"logfoldchanges", "pvals_adj"},
}
MISSING = {"", "na", "nan", "null", "none", "-", "n/a"}


def _number(value: str) -> float | None:
    v = value.strip()
    if v.lower() in MISSING:
        return None
    x = float(v)
    return None if math.isnan(x) else x


def read_table(text: str, filename: str = "") -> dict[str, Any]:
    """Parse a CSV or TSV differential table into standard columns. Raises ValueError with a readable message."""
    text = text.lstrip("﻿")
    if not text.strip():
        raise ValueError("the table is empty")
    first = text.split("\n", 1)[0]
    name = filename.lower()
    if name.endswith(".csv"):
        delimiter = ","
    elif name.endswith((".tsv", ".tab")):
        delimiter = "\t"
    else:
        delimiter = max(("\t", ",", ";"), key=first.count)
    rows = [r for r in csv.reader(io.StringIO(text), delimiter=delimiter) if any(c.strip() for c in r)]
    header, body = [h.strip() for h in rows[0]], rows[1:]
    if not body:
        raise ValueError("the table has a header but no rows")
    if len(body) > MAX_TABLE_ROWS:
        raise ValueError(f"the table has {len(body):,} rows; the demo takes up to {MAX_TABLE_ROWS:,}")
    # R's write.csv writes row names under an empty header one column short of the rows.
    if len(body[0]) == len(header) + 1:
        header = [""] + header

    lower = [h.lower() for h in header]
    mapping: dict[str, str] = {}
    for std, names in ALIASES.items():
        for alias in names:
            if alias in lower and header[lower.index(alias)] not in mapping.values():
                mapping[std] = header[lower.index(alias)]
                break
    notes = []
    if "feature" not in mapping:
        if header[0] in mapping.values():
            raise ValueError("no column of feature names found; name it 'feature' or 'gene', or put the names first")
        mapping["feature"] = header[0]
        notes.append(f"no feature-name column recognised; using the first column, '{header[0] or '(unnamed)'}'")
    for std in ("effect", "p"):
        if std not in mapping:
            raise ValueError(f"no column for '{std}' found. Recognised names: {', '.join(n for n in ALIASES[std] if n)}")

    index = {std: header.index(col) for std, col in mapping.items()}
    table: dict[str, list[Any]] = {std: [] for std in index}
    for i, row in enumerate(body, start=2):
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))
        for std, j in index.items():
            if std in ("feature", "symbol"):
                table[std].append(row[j].strip())
            else:
                try:
                    table[std].append(_number(row[j]))
                except ValueError:
                    raise ValueError(f"line {i}: '{row[j]}' in column '{header[j]}' is not a number") from None

    if "p_adj" not in table:
        idx = [i for i, p in enumerate(table["p"]) if p is not None]
        adjusted = enrich.bh([table["p"][i] for i in idx])
        table["p_adj"] = [None] * len(table["p"])
        for i, q in zip(idx, adjusted):
            table["p_adj"][i] = q
        notes.append("no adjusted p-value column; computed Benjamini-Hochberg from the p-values")

    columns = set(lower)
    tool = next((t for t, cols in TOOL_HINTS.items() if cols <= columns), None)
    return {"table": table, "mapping": mapping, "notes": notes, "n_rows": len(body), "tool_guess": tool,
            "unmapped": [h for h in header if h not in mapping.values()]}


def _split(text: str | None, separators: str = r"[,\n;]+") -> list[str]:
    return [s.strip() for s in re.split(separators, text or "") if s.strip()]


def _group_sizes(text: str | None) -> dict[str, int] | None:
    groups: dict[str, int] = {}
    for part in _split(text):
        m = re.fullmatch(r"(.+?)\s*[=:]\s*(\d+)", part)
        if not m:
            raise ValueError(f"group sizes: write each group as name=number, for example treated=6 (got '{part}')")
        groups[m.group(1).strip()] = int(m.group(2))
    return groups or None


def _crosstab(text: str | None) -> dict[str, dict[str, int]] | None:
    """A group-by-batch table: a header of batch names, then one row per group with its counts."""
    lines = [re.split(r"[\s,;]+", ln.strip()) for ln in (text or "").splitlines() if ln.strip()]
    if not lines:
        return None
    header, body = lines[0], lines[1:]
    if not body:
        raise ValueError("batch table: give a header of batch names, then one line per group with its counts")
    widths = {len(row) - 1 for row in body}
    if len(widths) > 1:
        raise ValueError("batch table: every group line needs the same number of counts")
    # The header is the batch names, or a label such as 'group' followed by them.
    batches = header if len(header) == widths.pop() else header[1:]
    out: dict[str, dict[str, int]] = {}
    for row in body:
        counts = row[1:]
        if len(batches) != len(counts):
            raise ValueError(f"batch table: the row for '{row[0]}' has {len(counts)} counts for {len(batches)} batches")
        try:
            out[row[0]] = {b: int(n) for b, n in zip(batches, counts)}
        except ValueError:
            raise ValueError(f"batch table: the counts for '{row[0]}' must be whole numbers") from None
    return out


def _int(value: Any, name: str, default: int | None = None) -> int | None:
    if value in (None, ""):
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a whole number") from None


def _float(value: Any, name: str, default: float | None = None) -> float | None:
    if value in (None, ""):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a number") from None


def build(form: dict[str, Any], table_text: str, filename: str = "", gmt_text: str = "",
          gmt_name: str = "") -> dict[str, Any]:
    """A one-contrast payload from an uploaded table and the form, with the same call written as Python."""
    parsed = read_table(table_text, filename)
    unit_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", (form.get("unit_id") or "contrast").strip()) or "contrast"
    definition = (form.get("definition") or "").strip()
    if not definition:
        raise ValueError("write the contrast definition in words, including the direction "
                         "(for example: treated - control; positive effect = higher in treated)")

    sample_level = form.get("observation_type", "sample") != "cell"
    group_sizes = _group_sizes(form.get("groups"))
    n_obs = _int(form.get("n_observations"), "observations tested")
    if n_obs is None:
        if not (sample_level and group_sizes):
            raise ValueError("give the number of observations tested, or the samples per group for a sample-level test")
        n_obs = sum(group_sizes.values())

    batch_mode = form.get("batch_in_design", "unknown")
    kwargs: dict[str, Any] = {
        "analysis_id": (form.get("analysis_id") or "web-demo").strip(),
        "question": (form.get("question") or "").strip() or None,
        "tool": (form.get("tool") or parsed["tool_guess"] or "unspecified").strip(),
        "steps": _split(form.get("steps"), r"[\n;]+") or ["not described"],
        "n_observations": n_obs,
        "observation_type": "sample" if sample_level else "cell",
        "group_sizes": group_sizes,
        "units_are_replicates": sample_level,
        "batch_crosstab": None if batch_mode == "single" else _crosstab(form.get("batch")),
        "batch_in_design": {"yes": True, "no": False}.get(batch_mode),
        "single_batch": batch_mode == "single",
        "mean_variance": (form.get("mean_variance") or "").strip() or None,
        "n_features": _int(form.get("n_features"), "features per direction", 15),
        "fdr": _float(form.get("fdr"), "FDR", 0.05),
        "min_effect": _float(form.get("min_effect"), "minimum effect", 0.0),
        "exclude": (form.get("exclude") or "").strip() or None,
        "panel": _split(form.get("panel")),
    }
    if kwargs["exclude"]:
        try:
            re.compile(kwargs["exclude"])
        except re.error as exc:
            raise ValueError(f"exclude pattern: {exc}") from None
    if group_sizes:
        kwargs["sizes"] = {unit_id: {"n_per_group": group_sizes}}
    if filename:
        kwargs["source_tables"] = {unit_id: filename}

    gene_sets = None
    if gmt_text.strip():
        with tempfile.NamedTemporaryFile("w", suffix=".gmt", delete=False, encoding="utf-8") as f:
            f.write(gmt_text)
        gene_sets = enrich.read_gmt(f.name)
        Path(f.name).unlink()
        if not gene_sets:
            raise ValueError("the gene-set file has no sets: each line should be name, description, then genes, tab-separated")
        kwargs["gene_sets"] = gene_sets
        kwargs["gene_sets_name"] = gmt_name or "uploaded.gmt"

    payload = core.dge_payload({unit_id: parsed["table"]}, {unit_id: definition}, **kwargs)
    text = json.dumps(core.jsonable(payload), indent=2, ensure_ascii=False)
    return {
        "payload_text": text,
        "inspect": inspect(text),
        "table": {k: parsed[k] for k in ("mapping", "notes", "n_rows", "tool_guess", "unmapped")},
        "input_chars": len(table_text),
        "python": _python_call(unit_id, definition, kwargs, parsed["mapping"]),
    }


def _python_call(unit_id: str, definition: str, kwargs: dict[str, Any], mapping: dict[str, str]) -> str:
    """The same call written against the installed package, for copying into a pipeline."""
    # pandas names an unnamed first column (R row names) 'Unnamed: 0'
    rename = {(col or "Unnamed: 0"): std for std, col in mapping.items() if col != std}
    source = kwargs.get("source_tables", {}).get(unit_id, "my_table.tsv")
    defaults = {"observation_type": "sample", "units_are_replicates": True, "n_features": 15, "fdr": 0.05,
                "min_effect": 0.0, "single_batch": False}
    shown = {k: v for k, v in kwargs.items()
             if k not in ("gene_sets", "source_tables") and v not in (None, [], {}) and defaults.get(k, object()) != v}
    if "gene_sets" in kwargs:
        shown["gene_sets"] = "__GMT__"
    lines = [
        "import pandas as pd",
        "import llomics",
        "from llomics import enrich",
        "",
        f"df = pd.read_csv({source!r}, sep=None, engine='python')",
    ]
    if rename:
        lines.append(f"df = df.rename(columns={rename!r})")
    lines += [
        f"table = {{c: df[c].tolist() for c in {sorted(mapping)!r}}}",
        "",
        "payload = llomics.dge_payload(",
        f"    {{{unit_id!r}: table}},",
        f"    {{{unit_id!r}: {definition!r}}},",
    ]
    for key, value in shown.items():
        if value == "__GMT__":
            lines.append(f"    gene_sets=enrich.read_gmt({kwargs.get('gene_sets_name', 'sets.gmt')!r}),")
        else:
            lines.append(f"    {key}={value!r},")
    lines += [")", f"llomics.write_payload(payload, 'outputs/{kwargs['analysis_id']}.llomics.json')"]
    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------
# Example table for the build form
# ---------------------------------------------------------------------------------------------

INTERFERON = ("ISG15 MX1 MX2 OAS1 OAS2 OAS3 OASL IFIT1 IFIT2 IFIT3 IFIT5 IFI6 IFI27 IFI44 IFI44L IFIH1 RIGI "
              "STAT1 STAT2 IRF7 IRF9 RSAD2 HERC5 HERC6 USP18 BST2 XAF1 EPSTI1 PARP9 PARP14 DTX3L GBP1 GBP4 "
              "EIF2AK2 SAMD9L CMPK2 PLSCR1 SP110 TRIM22 LY6E IFITM1 IFITM3 ZC3HAV1 RTP4 LGALS9 TNFSF10 SAMD9 "
              "DDX60 HELZ2 IFI35 PSMB8 PSMB9 TAP1 NMI SP100 ISG20").split()


def example() -> dict[str, Any]:
    """A simulated limma topTable for an interferon response, and form values that describe it.

    Simulated, not measured: about 50 interferon-stimulated genes are shifted up in 6 treated versus
    6 control samples; every other gene has no true effect. Deterministic, so the page and the tests
    see the same table.
    """
    rng = random.Random(20261119)
    background = (ROOT / "benchmark" / "data" / "gene_symbols.txt").read_text(encoding="utf-8").split()
    genes = list(dict.fromkeys(background + INTERFERON + ["ACTB", "GAPDH"]))
    responders = set(INTERFERON[:50])
    rows = []
    for g in genes:
        ave = min(14.0, max(2.0, rng.gauss(8.0, 1.8)))
        se = 0.10 + 0.25 * rng.random() + (0.15 if ave < 5 else 0.0)
        true = rng.uniform(1.2, 3.2) if g in responders else 0.0
        effect = true + rng.gauss(0.0, se)
        t = effect / se
        p = max(math.erfc(abs(t) / math.sqrt(2)), 1e-300)
        rows.append([g, effect, ave, t, p])
    adjusted = enrich.bh([r[4] for r in rows])
    rows = sorted((r + [q] for r, q in zip(rows, adjusted)), key=lambda r: r[4])
    out = io.StringIO()
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val"])
    for g, effect, ave, t, p, q in rows:
        w.writerow([g, f"{effect:.4g}", f"{ave:.4g}", f"{t:.4g}", f"{p:.4g}", f"{q:.4g}"])
    return {
        "filename": "example_interferon_toptable.tsv",
        "text": out.getvalue(),
        "form": {
            "analysis_id": "example-interferon",
            "question": "Which genes respond to interferon-beta?",
            "unit_id": "ifn_vs_control",
            "definition": "interferon-beta 6 h - untreated, primary fibroblasts, adjusted for batch; "
                          "positive effect = higher after interferon",
            "tool": "limma",
            "steps": "filterByExpr\nlmFit ~ 0 + group + batch\ncontrasts.fit ifn - control\neBayes(trend = TRUE)",
            "observation_type": "sample",
            "groups": "ifn=6, control=6",
            "n_observations": "",
            "batch": "        b1 b2\nifn     3  3\ncontrol 3  3",
            "batch_in_design": "yes",
            "mean_variance": "trend",
            "panel": "ISG15, MX1, ACTB, GAPDH",
            "n_features": 15,
            "fdr": 0.05,
            "min_effect": 0,
            "exclude": "",
        },
    }


PUBLIC: dict[str, Callable[..., Any]] = {
    "info": info,
    "gallery": gallery,
    "inspect": inspect,
    "build": build,
    "example": example,
}
