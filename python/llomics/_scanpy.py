"""Scanpy wrapper: run rank_genes_groups and build an llomics payload from the result."""

from __future__ import annotations

import json
import re
from importlib.metadata import version
from pathlib import Path
from typing import Any, Sequence

from ._core import (
    CONTRACT_VERSION, DEFAULT_EXCLUDE, P_FLOOR, apply_trust_overrides, check, jsonable, p_histogram,
    p_shape, sig, worst, write_payload,
)

QC_KEYS = ("n_genes_by_counts", "total_counts", "pct_counts_mt")
HARVEST = ("log1p", "hvg", "pca", "neighbors", "leiden", "louvain", "umap")


def rank_genes_groups(
    adata: Any,
    groupby: str,
    *,
    n_features: int = 15,
    panel: Sequence[str] = (),
    question: str | None = None,
    analysis_id: str | None = None,
    out: str | Path | None = None,
    fdr: float = 0.05,
    min_logfc: float = 0.0,
    exclude: str | None = DEFAULT_EXCLUDE,
    sample_key: str | None = None,
    batch_key: str | None = None,
    doublet_key: str | None = "predicted_doublet",
    qc_keys: Sequence[str] = QC_KEYS,
    steps: Sequence[str] | None = None,
    trust: dict[str, Any] | None = None,
    rerun: bool = True,
    key: str = "rank_genes_groups",
    write_table: bool = True,
    p_floor: float = P_FLOOR,
    **kwargs: Any,
) -> dict[str, Any]:
    """Run ``sc.tl.rank_genes_groups`` and return an llomics payload built from its result.

    Extra keyword arguments (``method``, ``reference``, ``groups``, ``use_raw``, ``layer`` ...) are
    passed to ``sc.tl.rank_genes_groups``. With ``reference="rest"`` (the default) each group is a
    cluster unit with its up-regulated markers. With ``reference="<group>"`` each comparison is a
    contrast unit with up and down features.

    Parameters
    ----------
    n_features
        Features per unit: top markers per cluster, or per direction for a contrast. Never padded.
    panel
        Fixed features reported for every unit whether or not they rank. A feature that is not
        in the matrix is written as null.
    question
        The scientific question, in a sentence. It is sent to the model with every unit.
    out
        Path for the payload JSON. The full marker table is written next to it.
    fdr, min_logfc, exclude
        Eligibility filter for the top features: adjusted p below ``fdr``, absolute log fold change
        above ``min_logfc``, and feature name not matching the ``exclude`` regex.
    sample_key, batch_key, doublet_key
        ``adata.obs`` columns holding the biological sample, the batch, and the predicted-doublet
        flag. They feed the trust checks; a check whose column is missing is recorded as not_run.
    trust
        Declarations for checks the wrapper cannot compute, for example
        ``{"cluster_stability": ("pass", "stable across resolutions 0.4-1.0")}``.
    rerun
        False reuses the result already stored in ``adata.uns[key]`` instead of running the test.
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc

    if rerun:
        sc.tl.rank_genes_groups(adata, groupby, pts=True, key_added=key, **kwargs)
    elif key not in adata.uns:
        raise KeyError(f"adata.uns['{key}'] not found: run with rerun=True")

    params = dict(adata.uns[key]["params"])
    if params.get("groupby") != groupby:
        raise ValueError(f"adata.uns['{key}'] was computed for groupby='{params.get('groupby')}', not '{groupby}'")
    method = params.get("method")
    if method == "logreg":
        raise NotImplementedError("method='logreg' returns no p-values; use wilcoxon or a t-test")
    reference = params.get("reference", "rest")
    use_raw = params.get("use_raw")
    if use_raw is None:
        use_raw = adata.raw is not None
    groups = list(adata.uns[key]["names"].dtype.names)
    obs = adata.obs
    labels = obs[groupby].astype(str)
    pattern = re.compile(exclude) if exclude else None

    def table(group: str) -> Any:
        df = sc.get.rank_genes_groups_df(adata, group=group, key=key)
        df = df[df["pvals"].notna()].copy()
        df["excluded"] = df["names"].str.contains(pattern) if pattern else False
        return df

    def rows(df: Any) -> list[dict[str, Any]]:
        out_rows = []
        for r in df.itertuples(index=False):
            row = {"feature": r.names, "effect": sig(r.logfoldchanges), "stat": sig(r.scores),
                   "p_adj": sig(max(r.pvals_adj, p_floor))}
            # scanpy reports the reference fraction only for reference="rest"
            if hasattr(r, "pct_nz_group"):
                row["pct_in"] = sig(r.pct_nz_group, 3)
            if hasattr(r, "pct_nz_reference"):
                row["pct_out"] = sig(r.pct_nz_reference, 3)
            out_rows.append(row)
        return out_rows

    tables = {g: table(g) for g in groups}
    n_floored = int(sum((t["pvals_adj"] < p_floor).sum() for t in tables.values()))

    filters = [
        {"rule": f"pvals_adj < {fdr} ({params.get('corr_method')})", "why": "FDR control"},
        {"rule": f"top {n_features} by test statistic; never padded", "why": "payload size"},
    ]
    if min_logfc > 0:
        filters.insert(1, {"rule": f"|logfoldchanges| > {min_logfc}", "why": "minimum effect size"})
    if pattern:
        filters.insert(0, {"rule": f"top features exclude names matching {exclude}",
                           "why": "technical signal; applied to the feature list only, not to the test"})

    trust_checks: list[dict[str, Any]]
    if reference == "rest":
        units, trust_checks, dataset_reference = _cluster_units(
            adata, labels, groups, tables, rows, panel, use_raw, params, n_features, fdr, min_logfc,
            sample_key, batch_key, doublet_key, qc_keys, np, pd, sc)
        payload_method = "clusters"
    else:
        units, trust_checks, dataset_reference = _contrast_units(
            adata, labels, groups, tables, rows, panel, reference, groupby, params, n_features, fdr,
            min_logfc, sample_key, batch_key)
        payload_method = "dge"

    source_tables: dict[str, str] = {"anndata": f"adata.uns['{key}']"}
    if out is not None and write_table:
        table_path = Path(str(out).removesuffix(".json") + ".table.tsv.gz")
        table_path.parent.mkdir(parents=True, exist_ok=True)
        pd.concat([t.assign(group=g) for g, t in tables.items()]).to_csv(table_path, sep="\t", index=False)
        source_tables["table"] = str(table_path)

    harvested = steps is None
    run_steps = list(steps) if steps is not None else _harvest_steps(adata, key)
    payload = {
        "llomics": CONTRACT_VERSION,
        "method": payload_method,
        "run": {
            "analysis_id": analysis_id or f"{groupby}-{method}-{reference}",
            "question": question,
            "tool": f"scanpy {version('scanpy')} sc.tl.rank_genes_groups(method='{method}', reference='{reference}', use_raw={use_raw})",
            "steps": run_steps,
            "steps_source": ("harvested from adata.uns; cell and gene QC filters are not recorded there"
                             if harvested else "supplied by the analyst"),
            "filters": filters,
            "p_floor": p_floor,
            "n_p_values_floored": n_floored,
            "effect_note": "scanpy log2 fold changes are computed from group means and become very large when a feature has no counts in one group; read them with pct_in and pct_out",
            "panel_units": ("mean expression across the unit's cells, in the units of the tested matrix"
                            if payload_method == "clusters" else "same columns as the feature rows"),
            "source_tables": source_tables,
        },
        "dataset": {
            "n_observations": int(adata.n_obs),
            "observation_type": "cell",
            "n_biological_samples": int(obs[sample_key].nunique()) if sample_key else None,
            "n_features": int(adata.raw.n_vars if use_raw else adata.n_vars),
            "reference": dataset_reference,
        },
        "trust": apply_trust_overrides(trust_checks, trust),
        "units": units,
    }
    payload = jsonable(payload)
    if out is not None:
        path = write_payload(payload, out)
        print(f"llomics: wrote {path} ({len(units)} units)")
    return payload


def _harvest_steps(adata: Any, key: str) -> list[str]:
    steps = []
    for name in HARVEST:
        if name in adata.uns:
            entry = adata.uns[name]
            params = entry.get("params", entry) if hasattr(entry, "get") else {}
            flat = {k: v for k, v in jsonable(dict(params)).items() if not isinstance(v, (dict, list))}
            steps.append(f"{name}: {json.dumps(flat)}" if flat else name)
    flat = {k: v for k, v in jsonable(dict(adata.uns[key]["params"])).items() if not isinstance(v, (dict, list))}
    steps.append(f"rank_genes_groups: {json.dumps(flat)}")
    return steps


def _eligible(df: Any, fdr: float, min_logfc: float) -> Any:
    return (df["pvals_adj"] < fdr) & (df["logfoldchanges"].abs() > min_logfc) & ~df["excluded"]


def _cluster_units(adata, labels, groups, tables, rows, panel, use_raw, params, n_features, fdr, min_logfc,
                   sample_key, batch_key, doublet_key, qc_keys, np, pd, sc):
    obs = adata.obs
    n_total = int(adata.n_obs)
    qc_present = [k for k in qc_keys if k in obs.columns]
    has_doublets = bool(doublet_key) and doublet_key in obs.columns
    doublet = obs[doublet_key].astype(bool) if has_doublets else None

    reference: dict[str, Any] = {f"median_{k}": sig(obs[k].median()) for k in qc_present}
    reference["definition"] = "medians and fractions across all cells in the dataset"
    if has_doublets:
        reference["pct_predicted_doublet"] = sig(100 * doublet.mean(), 3)

    var_names = adata.raw.var_names if use_raw else adata.var_names
    present = [g for g in panel if g in var_names]
    panel_means = None
    if present:
        expr = sc.get.obs_df(adata, keys=present, use_raw=use_raw, layer=params.get("layer"))
        panel_means = expr.groupby(labels.values, observed=True).mean()

    nearest: dict[str, str] = {}
    if "X_pca" in adata.obsm and len(groups) > 2:
        cent = pd.DataFrame(adata.obsm["X_pca"]).groupby(labels.values, observed=True).mean()
        corr = np.corrcoef(cent.values)
        np.fill_diagonal(corr, -np.inf)
        nearest = {str(g): str(cent.index[int(i)]) for g, i in zip(cent.index, corr.argmax(axis=1))}

    batches = obs[batch_key].astype(str) if batch_key else None
    batch_share = batches.value_counts(normalize=True) if batch_key else None
    batch_flags = []

    units = []
    for g in groups:
        mask = (labels == g).values
        df = tables[g]
        up = df[_eligible(df, fdr, min_logfc) & (df["logfoldchanges"] > 0)].sort_values("scores", ascending=False)
        top = up.head(n_features)

        facts: dict[str, Any] = {}
        for k in qc_present:
            med = obs.loc[mask, k].median()
            facts[f"median_{k}"] = sig(med)
            ref = reference[f"median_{k}"]
            facts[f"{k}_vs_reference"] = sig(med / ref, 3) if ref else None
        if has_doublets:
            pct = 100 * doublet[mask].mean()
            facts["pct_predicted_doublet"] = sig(pct, 3)
            ref = reference["pct_predicted_doublet"]
            facts["doublet_fraction_vs_reference"] = sig(pct / ref, 3) if ref else None
        if sample_key:
            facts["n_samples_contributing"] = int(obs.loc[mask, sample_key].nunique())
        if batch_key:
            share = batches[mask].value_counts(normalize=True)
            facts["dominant_batch"] = str(share.index[0])
            facts["dominant_batch_fraction"] = sig(share.iloc[0], 3)
            if share.iloc[0] >= 0.9 and batch_share[share.index[0]] < 0.7:
                batch_flags.append(f"{g} ({share.iloc[0]:.0%} {share.index[0]})")
        if g in nearest:
            facts["nearest_cluster"] = nearest[g]

        panel_state = {
            p: (sig(panel_means.loc[g, p]) if panel_means is not None and p in present else None) for p in panel
        }
        units.append({
            "unit_id": str(g),
            "unit_type": "cluster",
            "size": {"n_cells": int(mask.sum()), "pct_of_total": sig(100 * mask.sum() / n_total, 3)},
            "facts": facts,
            "top_features": rows(top),
            "panel": panel_state,
            "tail": {
                "n_shown": int(len(top)),
                "n_eligible": int(len(up)),
                "rule": f"pvals_adj < {fdr} and logFC > {min_logfc}, ranked by test statistic",
            },
        })

    method = params.get("method")
    samples_note = (f"{adata.obs[sample_key].nunique()} biological sample(s) in obs['{sample_key}']"
                    if sample_key else "sample_key not supplied, so the number of biological samples is unknown")
    checks = [
        check("pseudoreplication", "warn",
              f"{method} test with cells as replicates; {samples_note}. P-values rank markers and are not evidence of a population-level difference."),
        check("doublet_detection", "pass", f"obs['{doublet_key}'] present; {reference['pct_predicted_doublet']}% of cells predicted doublets")
        if has_doublets else check("doublet_detection", "not_run", "no doublet column found; pass doublet_key"),
        check("cluster_stability", "not_run", "cannot be computed from one clustering; declare it with trust="),
    ]
    if not batch_key:
        checks.append(check("batch_integration", "not_run", "batch_key not supplied; declare a single batch with trust="))
    elif batches.nunique() == 1:
        checks.append(check("batch_integration", "not_applicable", f"single batch in obs['{batch_key}']"))
    elif batch_flags:
        checks.append(check("batch_integration", "warn", "clusters dominated by one batch: " + ", ".join(batch_flags)))
    else:
        checks.append(check("batch_integration", "pass", "no cluster is 90% or more from a batch that holds under 70% of cells"))
    return units, checks, reference


def _contrast_units(adata, labels, groups, tables, rows, panel, reference, groupby, params, n_features, fdr,
                    min_logfc, sample_key, batch_key):
    obs = adata.obs
    method = params.get("method")
    units, shapes, shape_notes, rep_notes, rep_statuses, batch_statuses, batch_notes = [], [], [], [], [], [], []

    for g in groups:
        df = tables[g]
        ok = _eligible(df, fdr, min_logfc)
        up = df[ok & (df["logfoldchanges"] > 0)].sort_values("scores", ascending=False)
        down = df[ok & (df["logfoldchanges"] < 0)].sort_values("scores", ascending=True)
        # Features with a statistic of exactly 0 (typically no counts in either group) carry p = 1
        # and say nothing about the test, so they are left out of the histogram and counted apart.
        informative = df[df["scores"] != 0]
        h = p_histogram(informative["pvals"].tolist())
        status, note = p_shape(h)
        if status == "warn":
            note += ". In sparse single-cell data, lowly expressed features give discrete, conservative p-values, which produces the same shape"
        shapes.append(status)
        shape_notes.append(f"{g} vs {reference}: {note}")
        in_g, in_ref = (labels == g).values, (labels == reference).values
        near = (df["pvals_adj"] >= fdr) & (df["pvals_adj"] < 2 * fdr)
        q = informative["logfoldchanges"].quantile([0.01, 0.25, 0.5, 0.75, 0.99]).tolist()

        size: dict[str, Any] = {"n_tested": int(len(df)), "n_cells": {str(g): int(in_g.sum()), str(reference): int(in_ref.sum())}}
        if sample_key:
            n_s = {str(g): int(obs.loc[in_g, sample_key].nunique()), str(reference): int(obs.loc[in_ref, sample_key].nunique())}
            size["n_samples"] = n_s
            rep_statuses.append("fail" if min(n_s.values()) < 2 else "warn")
            rep_notes.append(f"{g} vs {reference}: {n_s}")
        if batch_key:
            ct = {str(k): obs.loc[m, batch_key].astype(str).value_counts().to_dict() for k, m in ((g, in_g), (reference, in_ref))}
            shared = set(ct[str(g)]) & set(ct[str(reference)])
            n_batches = len(set(ct[str(g)]) | set(ct[str(reference)]))
            if n_batches == 1:
                batch_statuses.append("not_applicable")
            elif not shared:
                batch_statuses.append("fail")
                batch_notes.append(f"{g} vs {reference}: no batch contains both groups")
            else:
                batch_statuses.append("warn")
                batch_notes.append(f"{g} vs {reference}: {ct}")

        lookup = df.set_index("names")
        panel_state = {}
        for p in panel:
            if p in lookup.index:
                r = lookup.loc[p]
                panel_state[p] = {"effect": sig(r["logfoldchanges"]), "stat": sig(r["scores"]),
                                  "p_adj": sig(max(r["pvals_adj"], P_FLOOR))}
            else:
                panel_state[p] = None

        top_up, top_down = up.head(n_features), down.head(n_features)
        units.append({
            "unit_id": f"{g}_vs_{reference}",
            "unit_type": "contrast",
            "definition": f"{g} - {reference} in obs['{groupby}']; positive logFC = higher in {g}",
            "size": size,
            "facts": {
                "fdr": fdr,
                "n_hits": int(len(up) + len(down)),
                "n_up": int(len(up)),
                "n_down": int(len(down)),
                "n_near_miss": int(near.sum()),
                "n_uninformative": int(len(df) - len(informative)),
                "p_histogram": h,
                "logfc_quantiles": dict(zip(("q01", "q25", "q50", "q75", "q99"), (sig(v, 3) for v in q))),
            },
            "top_features": {"up": rows(top_up), "down": rows(top_down)},
            "panel": panel_state,
            "tail": {
                "n_shown": int(len(top_up) + len(top_down)),
                "n_eligible": int(len(up) + len(down)),
                "rule": f"pvals_adj < {fdr} and |logFC| > {min_logfc}, ranked by |test statistic|, up to {n_features} per direction",
            },
            "trust": [check("pvalue_distribution", status, note)],
        })

    cells_note = "cells were the test units, so p-values do not reflect variation between samples; use pseudobulk for inference"
    if sample_key:
        replicates = check("replicates", worst(rep_statuses), f"samples per group: {'; '.join(rep_notes)}. {cells_note}")
    else:
        replicates = check("replicates", "not_run", f"sample_key not supplied. {cells_note}")

    if not batch_key:
        batch = check("batch_confounding", "not_run", "batch_key not supplied")
    elif all(s == "not_applicable" for s in batch_statuses):
        batch = check("batch_confounding", "not_applicable", f"single batch in obs['{batch_key}']")
    else:
        real = [s for s in batch_statuses if s != "not_applicable"]
        batch = check("batch_confounding", worst(real),
                      f"{method} has no covariates, so batch is not adjusted. Cells per batch: " + " | ".join(batch_notes))

    corr = params.get("corr_method")
    checks = [
        replicates,
        check("pvalue_distribution", worst(shapes), " | ".join(shape_notes)),
        batch,
        check("multiple_testing", "pass", f"{corr} within the one comparison reported") if len(groups) == 1
        else check("multiple_testing", "warn", f"{corr} within each of {len(groups)} comparisons; not adjusted across them"),
        check("mean_variance", "not_applicable", "rank-based test") if method == "wilcoxon"
        else check("mean_variance", "warn", f"{method} on the supplied matrix with no mean-variance modelling"),
    ]
    return units, checks, {"n_contrasts": len(groups), "reference_group": str(reference)}
