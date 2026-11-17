#!/usr/bin/env python3
"""Regenerate the test fixtures. Deterministic; standard library only.

    python tests/fixtures/make_fixtures.py

Writes the input tables and argument files. The expected payloads (expected_*.json) are written by
running the tests with LLOMICS_UPDATE_GOLDEN=1, and should be reviewed before they are committed.
"""

from __future__ import annotations

import json
import math
import random
from pathlib import Path

HERE = Path(__file__).parent


def bh(pvals: list[float | None]) -> list[float | None]:
    idx = sorted((i for i, p in enumerate(pvals) if p is not None), key=lambda i: pvals[i])
    m = len(idx)
    out: list[float | None] = [None] * len(pvals)
    prev = 1.0
    for rank in range(m, 0, -1):
        i = idx[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        out[i] = prev
    return out


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    def cell(v):
        if v is None:
            return "NA"
        if isinstance(v, bool):
            return "TRUE" if v else "FALSE"
        return repr(v) if isinstance(v, float) else str(v)

    lines = ["\t".join(header)] + ["\t".join(cell(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def dge_table(rng: random.Random, n: int, n_up: int, n_down: int, rename: dict[int, str] | None = None) -> list[list]:
    rows = []
    for i in range(n):
        effect = rng.gauss(0, 0.35)
        if i < n_up:
            effect += rng.uniform(0.8, 2.5)
        elif i < n_up + n_down:
            effect -= rng.uniform(0.8, 2.5)
        se = rng.uniform(0.15, 0.6)
        stat = effect / se
        p = math.erfc(abs(stat) / math.sqrt(2))
        symbol = f"SYM{i + 1:04d}"
        if i % 37 == 0:
            symbol = ""                      # no readable name
        if i in (3, n_up + 2):
            symbol = f"MT-GENE{i}"           # matched by the exclude pattern
        if rename and i in rename:
            symbol = rename[i]               # a real gene name, so the technical panels have something to find
        rows.append([f"F{i + 1:04d}", symbol, effect, rng.uniform(2, 12), stat, p])
    # Edge cases
    rows[0][4], rows[0][5] = 60.0, 0.0       # underflowed p-value
    rows[1][4] = rows[2][4] = 30.0           # tie in the statistic
    rows[1][5] = rows[2][5] = math.erfc(30.0 / math.sqrt(2))
    for i in range(n - 12, n - 6):           # uninformative: statistic exactly 0, p = 1
        rows[i][2], rows[i][4], rows[i][5] = 0.0, 0.0, 1.0
    rows[n - 1][5] = None                    # not tested
    p_adj = bh([r[5] for r in rows])
    for r, q in zip(rows, p_adj):
        r.append(q)
    rows[n_up - 1][6] = None                 # adjusted p missing while raw p is present
    return rows


def main() -> None:
    rng = random.Random(20261117)
    header = ["feature", "symbol", "effect", "ave_expr", "stat", "p", "p_adj"]
    write_tsv(HERE / "dge_A_vs_B.tsv", header,
              dge_table(rng, 400, 40, 25, {300: "FOS", 301: "Jun", 302: "EGR1", 303: "HSPA1A"}))
    write_tsv(HERE / "dge_C_vs_B.tsv", header,
              dge_table(rng, 400, 6, 3, {1: "XIST", 2: "RPS4Y1", 4: "Ddx3y", 5: "KDM5D", 6: "UTY", 7: "EIF1AY", 8: "TSIX"}))
    sets = {
        "SET_UP": [f"SYM{i:04d}" for i in range(5, 31)],
        "SET_DOWN": [f"sym{i:04d}" for i in range(41, 66)],     # lower case: matching ignores case
        "SET_RANDOM": [f"SYM{i:04d}" for i in range(200, 231)],
        "SET_SMALL": ["SYM0010", "SYM0011", "SYM0012"],          # below the minimum size: not tested
        "SET_ABSENT": [f"NOPE{i}" for i in range(20)],           # no member measured: not tested
        "SET_CLUSTER": [f"G{i:03d}" for i in range(1, 13)],
        "SET_OTHER": [f"G{i:03d}" for i in range(100, 121)],
    }
    (HERE / "sets.gmt").write_text(
        "".join(f"{name}\tfixture\t" + "\t".join(members) + "\n" for name, members in sets.items()),
        encoding="utf-8", newline="\n")
    (HERE / "dge_args.json").write_text(json.dumps({
        "definitions": {"A_vs_B": "A - B, adjusted for batch; positive effect = higher in A",
                        "C_vs_B": "C - B, adjusted for batch; positive effect = higher in C"},
        "analysis_id": "fixture-dge",
        "question": "Which features differ between the groups?",
        "tool": "fixture 1.0",
        "steps": ["simulated table"],
        "n_observations": 11,
        "group_sizes": {"A": 4, "B": 4, "C": 3},
        "batch_crosstab": {"A": {"b1": 2, "b2": 2}, "B": {"b1": 2, "b2": 2}, "C": {"b1": 3}},
        "batch_in_design": True,
        "mean_variance": "voom",
        "n_features": 5,
        "fdr": 0.05,
        "min_effect": 0.2,
        "exclude": "^MT-",
        "panel": ["F0002", "SYM0050", "F0400", "NOT_MEASURED"],
        "sizes": {"A_vs_B": {"n_per_group": {"A": 4, "B": 4}}, "C_vs_B": {"n_per_group": {"C": 3, "B": 4}}},
        "contrast_weights": {"A_vs_B": {"A": 1, "B": -1}, "C_vs_B": {"C": 1, "B": -1}},
        "filters": [{"rule": "expression filter", "why": "low counts"}],
        "source_tables": {"A_vs_B": "dge_A_vs_B.tsv", "C_vs_B": "dge_C_vs_B.tsv"},
        "trust": {"mean_variance": {"status": "pass", "note": "voom with quality weights"}},
    }, indent=2) + "\n", encoding="utf-8", newline="\n")

    # Clusters: three clusters, marker tables without a statistic for one variant and with for another
    rng = random.Random(7)
    clusters = {"c0": 120, "c1": 60, "c2": 20}
    marker_rows = []
    for ci, cl in enumerate(clusters):
        for i in range(150):
            effect = rng.gauss(0, 0.5) + (rng.uniform(1, 4) if i < 12 + 4 * ci else 0)
            stat = effect * rng.uniform(2, 6)
            p = math.erfc(abs(stat) / math.sqrt(2))
            name = f"G{i + 1:03d}" if i != 5 else f"RPL{i}"
            marker_rows.append([cl, name, effect, stat, p, rng.uniform(0.2, 1.0), rng.uniform(0.0, 0.5)])
    for cl in clusters:
        idx = [k for k, r in enumerate(marker_rows) if r[0] == cl]
        for k, q in zip(idx, bh([marker_rows[k][4] for k in idx])):
            marker_rows[k].append(q)
    write_tsv(HERE / "markers.tsv", ["cluster", "feature", "effect", "stat", "p", "pct_in", "pct_out", "p_adj"], marker_rows)

    obs_rows = []
    for cl, n in clusters.items():
        for j in range(n):
            batch = "b2" if cl == "c2" and j < 19 else ("b1" if j % 2 else "b3")
            obs_rows.append([cl, rng.randint(800, 4000) * (2 if cl == "c2" else 1), round(rng.uniform(2, 12), 3),
                             (j % 9 == 0) or (cl == "c2" and j % 2 == 0), batch, f"s{j % 3}"])
    write_tsv(HERE / "obs.tsv", ["label", "n_genes", "pct_mt", "doublet", "batch", "sample"], obs_rows)
    (HERE / "clusters_args.json").write_text(json.dumps({
        "analysis_id": "fixture-clusters",
        "question": "Which cell type is each cluster?",
        "tool": "fixture 1.0",
        "steps": ["simulated markers"],
        "test": "wilcoxon test",
        "n_features": 4,
        "fdr": 0.05,
        "min_effect": 0.25,
        "exclude": "^RP[LS]\\d",
        "panel_values": {cl: {"CD3D": round(0.5 + k, 3), "MS4A1": 0.0, "NOT_MEASURED": None}
                         for k, cl in enumerate(clusters)},
        "extra_facts": {"c0": {"nearest_cluster": "c1"}, "c1": {"nearest_cluster": "c0"}, "c2": {"nearest_cluster": "c0"}},
        "n_dataset_features": 150,
        "source_tables": {"table": "markers.tsv"},
    }, indent=2) + "\n", encoding="utf-8", newline="\n")
    print("fixtures written to", HERE)


if __name__ == "__main__":
    main()
