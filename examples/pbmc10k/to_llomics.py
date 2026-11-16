#!/usr/bin/env python3
"""Convert the pbmc10k-jev-annotation cluster summaries into an llomics payload (contract v0.1).

Usage:
    python to_llomics.py path/to/cluster_summaries.json payload.json

The source file is the output of 01_preprocess_pbmc.py in the pbmc10k-jev-annotation repository.
The run card below is transcribed from that script. Fields the source does not record are written
as null or not_run, never guessed.
"""

from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path
from typing import Any

P_FLOOR = 1e-300


def build_payload(clusters: list[dict[str, Any]]) -> dict[str, Any]:
    # Reference as used in the case study: median of the per-cluster medians.
    reference = {
        "median_genes": statistics.median(c["median_genes"] for c in clusters),
        "median_pct_mt": statistics.median(c["median_pct_mt"] for c in clusters),
        "median_pct_predicted_doublet": statistics.median(c["pct_predicted_doublet"] for c in clusters),
        "definition": "median of per-cluster medians (each cluster weighted equally)",
    }

    units = []
    n_floored = 0
    for c in clusters:
        rows = []
        for m in c["top_15_markers"]:
            p_adj = m["pval_adj"]
            if p_adj < P_FLOOR:
                p_adj = P_FLOOR
                n_floored += 1
            rows.append(
                {"feature": m["gene"], "effect": m["logfc"], "stat": m["wilcoxon_score"], "p_adj": p_adj}
            )
        units.append(
            {
                "unit_id": c["cluster_id"],
                "unit_type": "cluster",
                "size": {"n_cells": c["cell_count"], "pct_of_total": c["pct_of_total"]},
                "facts": {
                    "median_genes": c["median_genes"],
                    "median_umis": c["median_umis"],
                    "median_pct_mt": c["median_pct_mt"],
                    "median_doublet_score": c["median_doublet_score"],
                    "pct_predicted_doublet": c["pct_predicted_doublet"],
                    "genes_vs_reference": round(c["median_genes"] / reference["median_genes"], 2),
                    "pct_mt_minus_reference": round(c["median_pct_mt"] - reference["median_pct_mt"], 2),
                    "doublet_fraction_vs_reference": round(
                        c["pct_predicted_doublet"] / reference["median_pct_predicted_doublet"], 2
                    ),
                },
                "top_features": rows,
                "panel": c["canonical_marker_expression"],
                "tail": {
                    "n_shown": len(rows),
                    "n_eligible": None,  # the source keeps the top 15 and does not record how many passed
                    "rule": "adj. p < 0.05 and logFC > 0, ranked by adj. p then logFC, top 15",
                },
            }
        )

    return {
        "llomics": "0.1",
        "method": "clusters",
        "run": {
            "analysis_id": "pbmc10k-v3-leiden-r0.6",
            "question": "Which canonical PBMC population is each Leiden cluster, and can the label be trusted?",
            "tool": "scanpy sc.tl.rank_genes_groups(method='wilcoxon', use_raw=True), one cluster vs rest",
            "steps": [
                "10x Genomics pbmc_10k_v3 filtered feature-barcode matrix",
                "cells: n_genes >= 200, n_genes <= 5000, pct_counts_mt < 15; genes: min_cells = 3",
                "Scrublet on raw counts, expected_doublet_rate = 0.06, 30 PCs, exact neighbours",
                "normalize_total(1e4), log1p",
                "2000 highly variable genes, scale(max_value=10), PCA (arpack)",
                "neighbors(n_neighbors=15, n_pcs=30), Leiden(resolution=0.6), UMAP",
            ],
            "filters": [
                {"rule": "marker list excludes ^(MT-|RP[LS][0-9AP]|MALAT1$)",
                 "why": "technical signal; applied to the marker list only, not to clustering"},
                {"rule": "pvals_adj < 0.05 and logfoldchanges > 0", "why": "significant, up-regulated markers only"},
                {"rule": "top 15 by adjusted p then logFC; never padded", "why": "payload size"},
            ],
            "p_floor": P_FLOOR,
            "n_p_values_floored": n_floored,
            "panel_units": "mean log-normalised expression across the cluster's cells",
            "source_tables": {"anndata": "input/processed/pbmc10k_clustered.h5ad (uns['rank_genes_groups'])"},
        },
        "dataset": {
            "n_observations": sum(c["cell_count"] for c in clusters),
            "observation_type": "cell",
            "n_biological_samples": 1,
            "n_features": None,  # post-QC gene count is not recorded in the source summaries
            "n_units": len(clusters),
            "reference": reference,
        },
        "trust": [
            {"check": "pseudoreplication", "status": "warn",
             "note": "One donor; cells are the replicates. Wilcoxon p-values rank markers and are not evidence of a population-level difference."},
            {"check": "doublet_detection", "status": "pass",
             "note": "Scrublet on raw counts. The 6% expected doublet rate is an assumption."},
            {"check": "cluster_stability", "status": "not_run"},
            {"check": "batch_integration", "status": "not_applicable", "note": "Single sample, single batch."},
        ],
        "units": units,
    }


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    clusters = json.loads(src.read_text(encoding="utf-8"))
    payload = build_payload(clusters)
    dst.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    sizes = [len(json.dumps(u, separators=(",", ":"))) for u in payload["units"]]
    print(
        f"Wrote {dst}: {len(sizes)} units, {min(sizes)}-{max(sizes)} chars per unit (compact), "
        f"{payload['run']['n_p_values_floored']} p-values clamped to {P_FLOOR:g}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
