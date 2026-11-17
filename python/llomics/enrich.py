"""Technical-panel signatures and gene-set enrichment, computed in code. Standard library only.

Both answer questions a model should not be left to guess from fifteen gene names: is the hit list
dominated by a known artefact, and which gene sets are over-represented among all the hits. The R
package mirrors this module (r/llomics/R/enrich.R).
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Mapping, Sequence

from .contract import CONTRACT, DEFAULTS

PANELS: dict[str, Any] = json.loads(Path(__file__).with_name("llomics-panels.json").read_text(encoding="utf-8"))["panels"]


def _sig(x: float, digits: int) -> float:
    return float(f"{x:.{digits}g}")


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) when drawing n from N of which K are marked."""
    if k <= 0:
        return 1.0
    upper = min(K, n)
    if k > upper:
        return 0.0

    def logc(a: int, b: int) -> float:
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)

    denom = logc(N, n)
    terms = [logc(K, i) + logc(N - K, n - i) - denom for i in range(k, upper + 1) if n - i <= N - K]
    if not terms:
        return 0.0
    top = max(terms)
    return min(1.0, math.exp(top) * sum(math.exp(t - top) for t in terms))


def bh(pvals: Sequence[float]) -> list[float]:
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    out = [1.0] * len(pvals)
    prev = 1.0
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * len(order) / rank)
        out[i] = prev
    return out


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    """Read gene sets from a GMT file: one set per line, name, description, then members, tab-separated."""
    sets: dict[str, list[str]] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3 and parts[0]:
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def _in_panel(label_upper: str, genes: set[str], pattern: re.Pattern[str] | None) -> bool:
    return label_upper in genes or bool(pattern and pattern.search(label_upper))


def technical_signatures(tested_labels: Sequence[str], ranked_hit_labels: Sequence[str],
                         panels: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Score each technical panel against the hits.

    ``tested_labels`` are the readable names of every tested feature; ``ranked_hit_labels`` are the
    significant ones, best first. Returns the facts to merge into a unit. When no panel feature was
    tested at all, the flag count is None: the question could not be asked, which is not a zero.
    """
    panels = PANELS if panels is None else panels
    cfg = CONTRACT["technical_flag"]
    tested = [str(v).upper() for v in tested_labels]
    hits = [str(v).upper() for v in ranked_hit_labels]
    top = hits[:cfg["top_n"]]
    signatures: dict[str, Any] = {}
    flags: list[str] = []
    measured_total = 0
    for name, panel in panels.items():
        genes = {g.upper() for g in panel.get("genes", [])}
        pattern = re.compile(panel["pattern"], re.IGNORECASE) if panel.get("pattern") else None
        n_measured = sum(_in_panel(v, genes, pattern) for v in tested)
        measured_total += n_measured
        n_hits = sum(_in_panel(v, genes, pattern) for v in hits)
        if n_hits == 0:
            continue
        top_share = sum(_in_panel(v, genes, pattern) for v in top) / len(top)
        p = max(hypergeom_sf(n_hits, len(tested), n_measured, len(hits)), DEFAULTS["p_floor"])
        signatures[name] = {"n_measured": n_measured, "n_hits": n_hits, "top_share": _sig(top_share, 3), "p": _sig(p, 4)}
        if n_hits >= cfg["min_hits"] and p < cfg["max_p"] and top_share >= cfg["min_top_share"]:
            flags.append(name)
    return {
        "n_panel_features_measured": measured_total,
        "technical_signatures": signatures,
        "technical_flags": flags,
        "n_technical_flags": len(flags) if measured_total > 0 else None,
    }


def enrichment(tested_labels: Sequence[str], ranked_hit_labels: Sequence[str],
               gene_sets: Mapping[str, Sequence[str]]) -> tuple[list[dict[str, Any]], int]:
    """Over-representation of gene sets among the hits, against the tested features.

    One-sided hypergeometric test per set, Benjamini-Hochberg across the sets tested. Returns the
    significant sets (best first, capped) and the number of sets tested.
    """
    cfg = CONTRACT["enrichment"]
    universe = {str(v).upper() for v in tested_labels}
    hits: list[tuple[str, str]] = []
    seen: set[str] = set()
    for v in ranked_hit_labels:
        up = str(v).upper()
        if up not in seen:
            seen.add(up)
            hits.append((up, str(v)))

    tested_sets = []
    for name, members in gene_sets.items():
        inside = {str(g).upper() for g in members} & universe
        if len(inside) < cfg["min_set_size"]:
            continue
        overlap = [original for up, original in hits if up in inside]
        p = hypergeom_sf(len(overlap), len(universe), len(inside), len(hits))
        tested_sets.append((name, len(inside), overlap, p))

    adjusted = bh([t[3] for t in tested_sets])
    rows = []
    for (name, n_set, overlap, p), p_adj in zip(tested_sets, adjusted):
        if p_adj < cfg["fdr"] and len(overlap) >= cfg["min_overlap"]:
            rows.append((p, name, {
                "set": name,
                "n_set": n_set,
                "n_overlap": len(overlap),
                "expected": _sig(len(hits) * n_set / len(universe), 3),
                "p_adj": _sig(max(p_adj, DEFAULTS["p_floor"]), 4),
                "genes": overlap[:cfg["max_genes"]],
            }))
    rows.sort(key=lambda t: (t[0], t[1]))
    return [r[2] for r in rows[:cfg["max_sets"]]], len(tested_sets)
