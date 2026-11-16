---
name: llomics-lasso
description: >
  Build an llomics payload from a lasso or elastic-net model (glmnet, cv.glmnet, scikit-learn
  LassoCV, ElasticNetCV, LogisticRegressionCV with an L1 penalty): one state per model at the
  chosen penalty with cross-validated performance against a null model, the selected features with
  selection frequency, and their correlated stand-ins. Use when a sparse model picked a signature
  from thousands of features and the question is whether the signature is stable and what it
  means.
---

# llomics lasso

The fit object holds a full regularisation path: a coefficient for every feature at every penalty
value. The decision concerns one point on that path: which features were selected, how well the
model predicts, and whether the selection would survive a different split of the data. Read
`llomics-core` first.

## Unit

One model at its chosen λ. When several outcomes are modelled, one unit per outcome.

## Extract

| Field | glmnet | scikit-learn |
| :--- | :--- | :--- |
| chosen penalty | `cv$lambda.min`, `cv$lambda.1se` | `alpha_` (LassoCV), `C_` (LogisticRegressionCV) |
| CV curve | `cv$cvm`, `cv$cvsd`, `cv$nzero` | `mse_path_`, `alphas_`; `scores_` |
| coefficients | `coef(cv, s = "lambda.1se")` | `coef_` |
| mixing | `alpha` argument (1 = lasso) | `l1_ratio_` |
| null model | `cv$cvm` at the largest λ (intercept only) | `DummyRegressor` or `DummyClassifier` under the same CV |

Selection frequency is not in the fit object. Refit on resamples (bootstrap or repeated
subsampling at half the data) at the chosen λ and count how often each feature is selected.

## Unit state

```json
{
  "unit_id": "relapse_lasso",
  "unit_type": "model",
  "definition": "logistic lasso, relapse within 2 years ~ 18,402 standardised expression features",
  "size": {"n": 184, "n_events": 61, "n_features": 18402, "n_selected": 14},
  "facts": {
    "alpha": 1.0,
    "lambda_rule": "lambda.1se",
    "cv": "10-fold, repeated 20 times, stratified",
    "metric": "AUC",
    "cv_metric": 0.74,
    "cv_metric_sd": 0.04,
    "null_metric": 0.50,
    "n_selected_at_lambda_min": 37,
    "n_stable": 6,
    "stability_threshold": 0.8,
    "n_resamples": 200
  },
  "top_features": [
    {"feature": "ENSG00000141510", "symbol": "TP53", "effect": -0.41, "selection_freq": 0.97,
     "stand_ins": [{"symbol": "CDKN1A", "r": 0.83}]},
    {"feature": "ENSG00000148773", "symbol": "MKI67", "effect": 0.33, "selection_freq": 0.52,
     "stand_ins": [{"symbol": "TOP2A", "r": 0.91}, {"symbol": "BIRC5", "r": 0.88}]}
  ],
  "panel": {"ESR1": {"effect": 0.0, "selection_freq": 0.08}, "ERBB2": null},
  "tail": {"n_shown": 14, "n_eligible": 14, "rule": "all non-zero coefficients at lambda.1se, ranked by selection frequency"}
}
```

- **`effect`** is the coefficient on standardised features. It is shrunk toward zero by design and
  is not an effect size estimate. Say so in `definition` or a trust note.
- **`selection_freq`** is the share of resamples in which the feature was selected. It is the
  nearest thing the lasso has to a measure of confidence in a feature.
- **`stand_ins`** lists unselected features highly correlated with a selected one. The lasso keeps
  one of a correlated group, more or less arbitrarily. `MKI67` at frequency 0.52 with `TOP2A` at
  r = 0.91 means "a proliferation gene", not "MKI67".
- **`panel`** holds features expected to be predictive. A zero coefficient there is informative,
  and `selection_freq` separates "never selected" from "sometimes selected". `null` means the
  feature was not among the candidates.
- **`null_metric`** is the same metric from an intercept-only model under the same CV. An AUC of
  0.74 is read against 0.50; a deviance is unreadable without it.
- Send all selected features when there are 50 or fewer.

## Required trust checks

| Check | Pass when | Why it matters |
| :--- | :--- | :--- |
| `nested_cv` | the reported performance comes from data not used to choose λ (nested CV or a held-out set) | Performance at the λ that minimised CV error on the same folds is optimistic |
| `leakage` | every data-dependent step (filtering by variance or by association with the outcome, scaling, imputation, batch correction) was fitted inside each training fold | Filtering features on the full data before CV is the most common way to obtain a signature that predicts nothing |
| `selection_stability` | resampling was run and frequencies are reported | A single lasso fit is one draw. With correlated omics features the selected set changes with the fold seed |
| `null_comparison` | the null-model metric is reported, and ideally a permuted-outcome distribution | Shows the model beats chance by more than its own CV noise |

Also record: class balance and the metric's sensitivity to it, whether samples in different folds
are independent (repeated measures from one patient must stay in one fold), and that no p-values
are reported unless post-selection inference was run.

## Judgments and decisions

| Belongs to code | Belongs to the model |
| :--- | :--- |
| which features pass the stability threshold | Do the stable features and their stand-ins share a biological theme? |
| whether CV performance exceeds null by a margin | Is a selected feature a plausible driver, or a marker of a covariate (sex, batch, tissue composition)? |
| overlap with other signatures | Which correlated group does each selected feature most plausibly represent? |

## Pitfalls

- Presenting the selected set as "the important genes". It is one sparse solution among many of
  near-equal fit.
- Reporting coefficients from unstandardised features, which rank by measurement scale.
- Comparing `lambda.min` and `lambda.1se` selections without saying which was used for the
  reported performance.
- For elastic net, `alpha` below 1 keeps correlated groups together. State it, because it changes
  how the selected set should be read.
