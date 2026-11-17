"""Shared helpers for the llomics tests. Standard library only."""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "tests" / "fixtures"
sys.path.insert(0, str(ROOT / "python"))

UPDATE_GOLDEN = os.environ.get("LLOMICS_UPDATE_GOLDEN") == "1"
TEXT_COLUMNS = {"feature", "symbol", "cluster", "label", "batch", "sample"}


def read_tsv(path: Path) -> dict[str, list[Any]]:
    """Read a fixture table into columns. 'NA' is missing; TRUE/FALSE are booleans."""
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    header, body = rows[0], rows[1:]

    def parse(col: str, value: str) -> Any:
        if col in TEXT_COLUMNS:
            return value
        if value == "NA":
            return None
        if value in ("TRUE", "FALSE"):
            return value == "TRUE"
        return float(value)

    return {col: [parse(col, row[i]) for row in body] for i, col in enumerate(header)}


def split_by(table: dict[str, list[Any]], key: str) -> dict[str, dict[str, list[Any]]]:
    out: dict[str, dict[str, list[Any]]] = {}
    for i, group in enumerate(table[key]):
        target = out.setdefault(group, {c: [] for c in table if c != key})
        for c in target:
            target[c].append(table[c][i])
    return out


def build_dge() -> dict[str, Any]:
    from llomics import core, enrich

    args = json.loads((FIXTURES / "dge_args.json").read_text(encoding="utf-8"))
    tables = {uid: read_tsv(FIXTURES / f"dge_{uid}.tsv") for uid in args["definitions"]}
    definitions = args.pop("definitions")
    return core.dge_payload(tables, definitions, gene_sets=enrich.read_gmt(FIXTURES / "sets.gmt"),
                            gene_sets_name="sets.gmt", **args)


def build_clusters() -> dict[str, Any]:
    from llomics import core, enrich

    args = json.loads((FIXTURES / "clusters_args.json").read_text(encoding="utf-8"))
    tables = split_by(read_tsv(FIXTURES / "markers.tsv"), "cluster")
    obs = read_tsv(FIXTURES / "obs.tsv")
    return core.clusters_payload(
        tables, labels=obs["label"], qc={"n_genes": obs["n_genes"], "pct_mt": obs["pct_mt"]},
        doublet=obs["doublet"], batch=obs["batch"], sample=obs["sample"],
        gene_sets=enrich.read_gmt(FIXTURES / "sets.gmt"), gene_sets_name="sets.gmt", **args)


def differences(a: Any, b: Any, path: str = "payload", rel: float = 1.5e-3) -> list[str]:
    """Structural comparison. Strings, keys and integers must match exactly; floats within `rel`."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = [f"{path}: keys differ: only left {sorted(set(a) - set(b))}, only right {sorted(set(b) - set(a))}"] \
            if set(a) != set(b) else []
        for k in sorted(set(a) & set(b)):
            out += differences(a[k], b[k], f"{path}.{k}", rel)
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path}: lengths differ: {len(a)} vs {len(b)}"]
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += differences(x, y, f"{path}[{i}]", rel)
        return out
    num = (int, float)
    if isinstance(a, num) and isinstance(b, num) and not isinstance(a, bool) and not isinstance(b, bool):
        if isinstance(a, int) and isinstance(b, int):
            return [] if a == b else [f"{path}: {a} != {b}"]
        scale = max(abs(a), abs(b))
        return [] if abs(a - b) <= rel * scale else [f"{path}: {a} != {b}"]
    return [] if a == b else [f"{path}: {a!r} != {b!r}"]


def golden(name: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Return the committed expected payload, rewriting it first when LLOMICS_UPDATE_GOLDEN=1."""
    path = FIXTURES / name
    if UPDATE_GOLDEN:
        path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    return json.loads(path.read_text(encoding="utf-8"))
