---
name: llomics-clusters
description: >
  Build an llomics payload from clustering plus marker ranking (Scanpy rank_genes_groups, Seurat
  FindAllMarkers, or any one-vs-rest marker table): one state per cluster with size, QC facts
  compared to the dataset, top markers, a canonical marker panel, and trust checks. Use for
  single-cell or spatial cluster annotation, cluster QC triage, and doublet or low-quality cluster
  review.
---

# llomics clusters

The full result is a cells × genes matrix plus a clusters × genes marker table. The decision is per
cluster: what is it, and can the label be trusted. Read `llomics-core` first.

This skill generalises the worked example in `examples/pbmc10k` (10,639 cells, 19 Leiden clusters,
a 51 KB payload from a 415 MB AnnData).

## Unit

One cluster. Send `dataset.reference` with every cluster so QC is judged against this dataset and
not against fixed cutoffs.

## Extract

| Field | Scanpy | Seurat |
| :--- | :--- | :--- |
| markers | `sc.tl.rank_genes_groups(adata, groupby, method="wilcoxon", pts=True)` then `sc.get.rank_genes_groups_df(adata, group)` | `FindAllMarkers(obj, only.pos = TRUE)` |
| effect | `logfoldchanges` | `avg_log2FC` |
| stat | `scores` | not returned; rank by `p_val_adj` then `avg_log2FC` |
| fraction expressing | `pct_nz_group`, `pct_nz_reference` | `pct.1`, `pct.2` |
| QC per cell | `adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]]` | `obj$nFeature_RNA`, `obj$nCount_RNA`, `obj$percent.mt` |
| doublets | Scrublet `doublet_score`, `predicted_doublet` | scDblFinder `scDblFinder.score`, `.class` |

## Wrappers

The llomics packages implement this skill:

```python
payload = llomics.rank_genes_groups(adata, "leiden", method="wilcoxon", n_features=15,
                                    panel=["CD3D", "MS4A1"], question="Which cell type is each cluster?",
                                    sample_key="donor", batch_key="run", out="markers.llomics.json")
```

```r
payload <- llomics_FindAllMarkers(obj, only.pos = TRUE, n_features = 15, panel = c("CD3D", "MS4A1"),
                                  question = "Which cell type is each cluster?",
                                  sample_col = "donor", batch_col = "run", out = "markers.llomics.json")
```

Both run the test, compute the facts and trust checks below, score the technical panels against
each cluster's markers, and add enrichment when `gene_sets` is given. Seurat reports no test
statistic, so its markers are ranked by p-value and then fold change, and its adjusted p-value is
a Bonferroni correction over all features.

## Unit state

```json
{
  "unit_id": "9",
  "unit_type": "cluster",
  "size": {"n_cells": 82, "pct_of_total": 0.77},
  "facts": {
    "median_genes": 3364.5,
    "median_umis": 12165.0,
    "median_pct_mt": 6.0,
    "pct_predicted_doublet": 14.63,
    "genes_vs_reference": 1.66,
    "pct_mt_minus_reference": -2.59,
    "doublet_fraction_vs_reference": 4.06,
    "nearest_cluster": "1"
  },
  "top_features": [
    {"feature": "PLD4", "effect": 7.028, "stat": 15.6, "p_adj": 3.626e-51, "pct_in": 0.95, "pct_out": 0.04}
  ],
  "panel": {"CD3D": 0.0695, "MS4A1": 0.1389, "LYZ": 0.8058, "CST3": 2.7975, "LILRA4": null},
  "tail": {"n_shown": 15, "n_eligible": 1840, "rule": "adj. p < 0.05 and logFC > 0, ranked by adj. p then logFC"}
}
```

The values are cluster 9 of the case study, shortened to one marker row. Four fields are
illustrative because the case study does not emit them yet: `pct_in`, `pct_out`,
`nearest_cluster`, and `tail.n_eligible`. `LILRA4` is shown as `null` to illustrate a panel gene
that was not measured.

- **`pct_in` and `pct_out`** separate a marker expressed in most cells of the cluster from one
  driven by a few high cells. Log fold change alone cannot.
- **`panel`** holds mean log-normalised expression of a fixed lineage panel. It shows the absence
  of other lineages, which is what distinguishes a clean cluster from a doublet cluster. A gene
  missing from the filtered matrix is `null`.
- **`nearest_cluster`** (by centroid correlation or PAGA connectivity) tells the reader what this
  cluster is most likely a subtype or a doublet of.
- **`facts` comparisons** are computed in code against `dataset.reference`.

## Required trust checks

| Check | Pass when | Why it matters |
| :--- | :--- | :--- |
| `pseudoreplication` | markers were tested across biological replicates (pseudobulk) | Cell-level tests treat cells from one donor as independent. With thousands of cells every expressed gene is "significant". Otherwise `warn`: p-values rank, they do not prove |
| `doublet_detection` | a doublet caller was run on raw counts; record the expected rate | The expected rate is an assumption and shifts every cluster's fraction |
| `cluster_stability` | clusters persist across resolutions or seeds | A cluster that appears at one resolution only may be an artefact of the parameter |
| `batch_integration` | batches are mixed within clusters, or there is one batch (`not_applicable`) | A cluster made of one batch is a batch effect until shown otherwise |

## Declare in `run.filters`

Marker exclusion patterns (mitochondrial, ribosomal, `MALAT1`), the significance and direction
cutoffs, and the cell and gene QC thresholds. State that excluded genes were removed from the
marker list only, or from clustering too.

## Judgments and decisions

| Belongs to code | Belongs to the model |
| :--- | :--- |
| doublet fraction against the reference | Which population do these markers define? |
| mitochondrial fraction and gene count against the reference | Do the markers form one coherent lineage program? |
| cluster size below a minimum | Are markers of mutually exclusive lineages co-expressed? |
| trust check statuses | Are the markers dominated by stress or cell-cycle genes? |

See `llomics-judge` for the gating pattern and for the case-study evidence on why the doublet
fraction belongs in code.

## Pitfalls

- With large clusters, adjusted p-values underflow to 0. In the case study 71 of 285 marker rows
  carried `pval_adj: 0.0`. Clamp to `run.p_floor` and rank by the statistic.
- A reference computed as the median of cluster medians weights a 26-cell cluster the same as a
  1,896-cell cluster. Say which reference is used; the median across cells is usually the one
  meant.
- Housekeeping genes with tiny fold changes rank high by p-value in large clusters. Rank by
  statistic with a minimum effect, or the top 15 fills with `EEF1A1` and `TPT1`.
