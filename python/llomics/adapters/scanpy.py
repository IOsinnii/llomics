"""Scanpy adapter: run rank_genes_groups, map its columns onto the standard table, call the core."""

from __future__ import annotations

import json
from importlib.metadata import version
from pathlib import Path
from typing import Any, Sequence

from .. import core, enrich
from ..contract import DEFAULTS

QC_KEYS = ("n_genes_by_counts", "total_counts", "pct_counts_mt")
HARVEST = ("log1p", "hvg", "pca", "neighbors", "leiden", "louvain", "umap")
COLUMNS = {"names": "feature", "logfoldchanges": "effect", "scores": "stat", "pvals": "p",
           "pvals_adj": "p_adj", "pct_nz_group": "pct_in", "pct_nz_reference": "pct_out"}
EFFECT_NOTE = ("scanpy log2 fold changes are computed from group means and become very large when a feature "
               "has no counts in one group; read them with pct_in and pct_out")


def rank_genes_groups(
    adata: Any,
    groupby: str,
    *,
    n_features: int | None = None,
    panel: Sequence[str] = (),
    question: str | None = None,
    analysis_id: str | None = None,
    out: str | Path | None = None,
    fdr: float | None = None,
    min_logfc: float = 0.0,
    exclude: str | None = DEFAULTS["exclude"],
    sample_key: str | None = None,
    batch_key: str | None = None,
    doublet_key: str | None = "predicted_doublet",
    qc_keys: Sequence[str] = QC_KEYS,
    steps: Sequence[str] | None = None,
    trust: dict[str, Any] | None = None,
    rerun: bool = True,
    key: str = "rank_genes_groups",
    write_table: bool = True,
    gene_sets: dict[str, Sequence[str]] | str | Path | None = None,
    gene_sets_name: str | None = None,
    technical_panels: dict[str, Any] | bool = True,
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
    gene_sets
        Gene sets for over-representation of the hits: a mapping of name to members, or the path
        to a GMT file.
    technical_panels
        True scores the shipped artefact panels (sex, stress, cell cycle ...) against the hits;
        False skips it; a mapping replaces the panels.
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc

    if isinstance(gene_sets, (str, Path)):
        gene_sets_name = gene_sets_name or Path(gene_sets).name
        gene_sets = enrich.read_gmt(gene_sets)
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
    groups = [str(g) for g in adata.uns[key]["names"].dtype.names]
    obs = adata.obs
    labels = obs[groupby].astype(str)
    if sample_key and sample_key not in obs.columns:
        raise KeyError(f"sample_key '{sample_key}' is not in adata.obs")
    if batch_key and batch_key not in obs.columns:
        raise KeyError(f"batch_key '{batch_key}' is not in adata.obs")

    frames = {g: sc.get.rank_genes_groups_df(adata, group=g, key=key) for g in groups}
    tables = {g: {COLUMNS[c]: df[c].tolist() for c in df.columns if c in COLUMNS} for g, df in frames.items()}

    source_tables: dict[str, str] = {"anndata": f"adata.uns['{key}']"}
    if out is not None and write_table:
        table_path = Path(str(out).removesuffix(".json") + ".table.tsv.gz")
        table_path.parent.mkdir(parents=True, exist_ok=True)
        pd.concat([df.assign(group=g) for g, df in frames.items()]).to_csv(table_path, sep="\t", index=False)
        source_tables["table"] = str(table_path)

    tool = (f"scanpy {version('scanpy')} sc.tl.rank_genes_groups(method='{method}', reference='{reference}', "
            f"use_raw={use_raw})")
    run_extra = {
        "steps_source": ("supplied by the analyst" if steps is not None else
                         "harvested from adata.uns; cell and gene QC filters are not recorded there"),
        "effect_note": EFFECT_NOTE,
    }
    run_steps = list(steps) if steps is not None else _harvest_steps(adata, key)
    common = dict(analysis_id=analysis_id or f"{groupby}-{method}-{reference}", question=question, tool=tool,
                  steps=run_steps, n_features=n_features, fdr=fdr, min_effect=min_logfc, exclude=exclude,
                  source_tables=source_tables, run_extra=run_extra, trust=trust, gene_sets=gene_sets,
                  gene_sets_name=gene_sets_name, technical_panels=technical_panels)

    if reference == "rest":
        qc = {k: obs[k].tolist() for k in qc_keys if k in obs.columns}
        has_doublets = bool(doublet_key) and doublet_key in obs.columns
        var_names = adata.raw.var_names if use_raw else adata.var_names
        present = [g for g in panel if g in var_names]
        panel_values: dict[str, dict[str, Any]] = {g: {p: None for p in panel} for g in groups}
        if present:
            expr = sc.get.obs_df(adata, keys=present, use_raw=use_raw, layer=params.get("layer"))
            means = expr.groupby(labels.values, observed=True).mean()
            for g in groups:
                for p in present:
                    panel_values[g][p] = core.sig(means.loc[g, p])

        extra_facts: dict[str, dict[str, Any]] = {}
        if "X_pca" in adata.obsm and len(groups) > 2:
            cent = pd.DataFrame(adata.obsm["X_pca"]).groupby(labels.values, observed=True).mean()
            corr = np.corrcoef(cent.values)
            np.fill_diagonal(corr, -np.inf)
            extra_facts = {str(g): {"nearest_cluster": str(cent.index[int(i)])}
                           for g, i in zip(cent.index, corr.argmax(axis=1))}

        payload = core.clusters_payload(
            tables, test=f"{method} test", labels=labels.tolist(), qc=qc,
            doublet=obs[doublet_key].astype(bool).tolist() if has_doublets else None,
            batch=obs[batch_key].astype(str).tolist() if batch_key else None,
            sample=obs[sample_key].astype(str).tolist() if sample_key else None,
            panel_values=panel_values, extra_facts=extra_facts,
            n_dataset_features=int(adata.raw.n_vars if use_raw else adata.n_vars), **common)
    else:
        reference = str(reference)
        in_ref = (labels == reference).values
        sizes, definitions, named_tables = {}, {}, {}
        sample_sizes: dict[str, int] = {}
        crosstab: dict[str, dict[str, int]] = {}
        for g in groups:
            uid = f"{g}_vs_{reference}"
            in_g = (labels == g).values
            named_tables[uid] = tables[g]
            definitions[uid] = f"{g} - {reference} in obs['{groupby}']; positive effect = higher in {g}"
            sizes[uid] = {"n_cells": {g: int(in_g.sum()), reference: int(in_ref.sum())}}
            for name, mask in ((g, in_g), (reference, in_ref)):
                if sample_key:
                    sample_sizes[name] = int(obs.loc[mask, sample_key].nunique())
                if batch_key:
                    crosstab[name] = {str(k): int(v) for k, v in obs.loc[mask, batch_key].astype(str).value_counts().items()}
            if sample_key:
                sizes[uid]["n_samples"] = {g: sample_sizes[g], reference: sample_sizes[reference]}
        payload = core.dge_payload(
            named_tables, definitions, n_observations=int(adata.n_obs), observation_type="cell",
            group_sizes=sample_sizes or None, units_are_replicates=False,
            batch_crosstab=crosstab or None, batch_adjustable=False,
            mean_variance="not_applicable" if method == "wilcoxon" else "none",
            correction=str(params.get("corr_method")), sizes=sizes, panel=panel,
            dataset_reference={"reference_group": reference}, **common)

    if out is not None:
        path = core.write_payload(payload, out)
        print(f"llomics: wrote {path} ({len(payload['units'])} units)")
    return payload


def _harvest_steps(adata: Any, key: str) -> list[str]:
    def flat(params: Any) -> dict[str, Any]:
        return {k: v for k, v in core.jsonable(dict(params)).items() if not isinstance(v, (dict, list))}

    steps = []
    for name in HARVEST:
        if name in adata.uns:
            entry = adata.uns[name]
            params = flat(entry.get("params", entry)) if hasattr(entry, "get") else {}
            steps.append(f"{name}: {json.dumps(params)}" if params else name)
    steps.append(f"rank_genes_groups: {json.dumps(flat(adata.uns[key]['params']))}")
    return steps
