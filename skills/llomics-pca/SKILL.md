---
name: llomics-pca
description: >
  Build an llomics payload from a PCA (prcomp, scikit-learn PCA, Scanpy sc.tl.pca): one state per
  component with variance explained, top positive and negative loadings, association with sample
  metadata, group positions, and outliers. Use when asking what each principal component
  represents, whether a component is biological or technical, or whether batch drives the main
  axes of variation.
---

# llomics PCA

The full result is two matrices: scores (samples × components) and loadings (features ×
components). Neither is readable. What a person wants to know about each component is what it
separates and which features drive it. Read `llomics-core` first.

## Unit

One component. Include components until cumulative variance passes a stated level or the scree
flattens, with a cap (10 is usually enough). Record the rule in `trust` under
`component_selection`.

## Extract

| Field | R `prcomp` | scikit-learn `PCA` | Scanpy |
| :--- | :--- | :--- | :--- |
| variance ratio | `fit$sdev^2 / sum(fit$sdev^2)` | `explained_variance_ratio_` | `adata.uns["pca"]["variance_ratio"]` |
| loadings | `fit$rotation` (features × PCs) | `components_.T` | `adata.varm["PCs"]` |
| scores | `fit$x` | `transform(X)` | `adata.obsm["X_pca"]` |

Check the denominator of the variance ratio. scikit-learn and Scanpy report it against total
variance even when only some components are computed. Truncated R solvers such as
`irlba::prcomp_irlba` return `sdev` for the kept components only; divide by `totalvar`, not by
`sum(sdev^2)`.

## Unit state

```json
{
  "unit_id": "PC2",
  "unit_type": "component",
  "facts": {
    "variance_explained": 0.143,
    "cumulative_variance": 0.521,
    "loading_concentration": 0.31,
    "orientation": "positive scores = batch B",
    "associations": [
      {"covariate": "batch", "kind": "categorical", "r2": 0.81, "p": 2e-9},
      {"covariate": "library_size", "kind": "continuous", "r2": 0.12, "p": 0.09},
      {"covariate": "treatment", "kind": "categorical", "r2": 0.04, "p": 0.41}
    ],
    "group_means": {"batch": {"A": -9.8, "B": 10.4}},
    "outliers": [{"sample": "S17", "robust_z": 4.6}]
  },
  "top_features": {
    "positive": [{"feature": "ENSG00000229807", "symbol": "XIST", "effect": 0.081}],
    "negative": [{"feature": "ENSG00000129824", "symbol": "RPS4Y1", "effect": -0.077}]
  },
  "panel": {"XIST": {"effect": 0.081, "rank": 1}, "MT-CO1": {"effect": 0.002, "rank": 9120}, "HBB": null},
  "tail": {"n_shown": 30, "n_eligible": 2000, "rule": "15 largest loadings per sign"}
}
```

Compute in code and store in `facts`:

- **`associations`**: for each metadata column, the share of the component's score variance it
  explains. Use R² from a one-way ANOVA for categorical covariates and squared Pearson r for
  continuous ones. Include technical covariates (batch, run date, library size, RIN, percent
  mitochondrial) next to the biological ones, sorted by R². This table is the most useful part of
  the payload.
- **`orientation`**: the sign of a component is arbitrary. State in words which samples sit on the
  positive side, so "positive loadings" has a meaning.
- **`loading_concentration`**: sum of squared loadings of the features shown. Near 1 means a few
  features define the component. Near 0 means the top list is a small sample of a diffuse axis and
  should not be read as "the genes of PC2".
- **`outliers`**: samples beyond a robust z of 3 on this component (median and MAD).
- **`panel`**: loading and rank of known technical and nuisance features: sex genes (`XIST`,
  `RPS4Y1`), mitochondrial and ribosomal genes, haemoglobin, cell-cycle genes.

## Required trust checks

| Check | Pass when | Why it matters |
| :--- | :--- | :--- |
| `scaling_declared` | transform, centring, scaling, and feature selection are stated in `run.steps` | Unscaled PCA follows the highest-variance features; scaled PCA weights noise features equally |
| `technical_covariate_association` | no retained component has high R² with a technical covariate; `warn` with the covariate named if one does | A component that is batch is not biology, and correlated biology cannot be separated from it |
| `outlier_influence` | no single sample contributes a dominant share of a component's score variance | One aberrant sample can create its own component |
| `component_selection` | the rule for how many components are reported is stated | Otherwise the reader cannot tell whether PC6 was dropped or never looked at |

## Judgments and decisions

| Belongs to code | Belongs to the model |
| :--- | :--- |
| which covariate has the highest R² | Do the positive loadings form a known program (cell cycle, interferon, hypoxia, sex)? |
| whether R² with batch exceeds a threshold | Do the positive and negative sides describe opposing states of one process? |
| which samples are outliers | Is this feature set more consistent with a technical source than with biology? |

## Pitfalls

- Reading loadings without the orientation statement.
- Interpreting the top 15 loadings of a diffuse component.
- A first component above half the variance in count data is often library size or a failed
  normalisation. Check its association with library size before interpreting it.
- UMAP and t-SNE coordinates are not components. They have no loadings and distances between
  groups are not interpretable; do not build this payload from them.
