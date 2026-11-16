---
name: llomics-plsda
description: >
  Build an llomics payload from PLS-DA, sPLS-DA, or OPLS-DA (mixOmics, ropls, scikit-learn
  PLSRegression on dummy-coded classes): a model state with cross-validated performance and a
  permutation test, and one state per component with high-VIP features, their direction, and a
  univariate cross-check. Use for metabolomics, proteomics, or lipidomics class discrimination
  when the question is which features separate the groups and whether the separation is real.
---

# llomics PLS-DA

PLS-DA will separate any two groups in a scores plot when features outnumber samples, including
groups assigned at random. The scores plot and the VIP list are therefore the least informative
parts of the result on their own. What decides whether they mean anything is cross-validated
performance and a permutation test. Read `llomics-core` first.

## Units

1. **The model**: classes, components, performance, permutation result. This unit carries the
   trust decision.
2. **Each component**: what it separates and which features drive it.

Interpret component units only if the model unit passes its gates.

## Extract

| Field | mixOmics | ropls | scikit-learn |
| :--- | :--- | :--- | :--- |
| fit | `plsda(X, Y, ncomp)` | `opls(X, y, predI, orthoI)` | `PLSRegression(n_components).fit(X, Y_dummy)` |
| CV performance | `perf(fit, validation = "Mfold", nrepeat)` → `$error.rate` (use `BER`) | `getSummaryDF(fit)`: `R2X(cum)`, `R2Y(cum)`, `Q2(cum)` | `cross_val_predict`, then balanced accuracy or AUC |
| permutation | not built in for PLS-DA; permute `Y` and rerun `perf` | `opls(..., permI = 1000)`, then `getSummaryDF(fit)`: `pR2Y`, `pQ2` | `permutation_test_score` |
| VIP | `vip(fit)` | `getVipVn(fit)` | compute from weights (below) |
| weights, loadings | `fit$loadings$X` | `getWeightMN`, `getLoadingMN` | `x_weights_`, `x_loadings_` |
| scores | `fit$variates$X` | `getScoreMN` | `x_scores_` |

VIP for feature j over A components, with p features, weights w, and SSY_a the Y sum of squares
explained by component a:

```text
VIP_j = sqrt( p * sum_a( SSY_a * (w_ja / ||w_a||)^2 ) / sum_a( SSY_a ) )
```

## Model unit state

```json
{
  "unit_id": "model",
  "unit_type": "model",
  "definition": "PLS-DA, case vs control, 312 metabolites, log-transformed and Pareto-scaled",
  "size": {"n_case": 28, "n_control": 31, "n_features": 312, "n_components": 2},
  "facts": {
    "r2x_cum": 0.34,
    "r2y_cum": 0.81,
    "q2_cum": 0.42,
    "r2y_minus_q2": 0.39,
    "cv": "7-fold, repeated 50 times, stratified",
    "cv_balanced_error": 0.21,
    "cv_auc": 0.86,
    "cv_confusion": {"case": {"case": 22, "control": 6}, "control": {"case": 6, "control": 25}},
    "permutations": 1000,
    "p_q2": 0.004,
    "p_accuracy": 0.006,
    "n_vip_above_1": 74
  },
  "top_features": [],
  "panel": {},
  "tail": {"n_shown": 0, "n_eligible": 0, "rule": "features are reported per component"}
}
```

## Component unit state

```json
{
  "unit_id": "comp1",
  "unit_type": "component",
  "facts": {
    "r2x": 0.22,
    "r2y": 0.63,
    "q2": 0.38,
    "orientation": "positive scores = case",
    "class_means": {"case": 2.9, "control": -2.6},
    "vip_quantiles": {"q50": 0.71, "q90": 1.38, "q99": 2.2}
  },
  "top_features": [
    {"feature": "HMDB0000161", "name": "L-Alanine", "vip": 2.41, "effect": 0.19,
     "higher_in": "case", "log2fc": 0.62, "p_adj": 0.003, "vip_cv_sd": 0.21}
  ],
  "panel": {"Glucose": {"vip": 0.44, "higher_in": "control", "p_adj": 0.61}, "Lactate": null},
  "tail": {"n_shown": 20, "n_eligible": 74, "rule": "VIP > 1, ranked by VIP"}
}
```

- **`higher_in`** is computed from group means in code. Loading signs depend on the arbitrary
  orientation of the component and are misread constantly. `effect` holds the weight for those who
  need it.
- **Univariate cross-check.** `log2fc` and `p_adj` from a simple per-feature test sit next to the
  VIP. A high-VIP feature with no univariate difference is contributing through correlation
  structure, which is legitimate but should be read with more care.
- **`vip_cv_sd`**: spread of the VIP across CV folds or resamples. A VIP of 1.2 ± 0.5 is not a
  finding.
- **`cv_confusion`** comes from out-of-fold predictions only.
- For OPLS-DA, report the predictive component as above and the orthogonal components as separate
  units with their metadata associations (see `llomics-pca`): they often hold batch or run order.

## Required trust checks

| Check | Pass when | Why it matters |
| :--- | :--- | :--- |
| `permutation_test` | class labels were permuted at least several hundred times, rerunning the full CV, and the observed Q² or accuracy sits in the tail | The only direct evidence that the separation is not what PLS-DA produces from noise. `not_run` should block interpretation |
| `cv_performance` | error, balanced error, or AUC is reported from held-out predictions | Fitted R²Y and the scores plot describe the training data only |
| `overfit_gap` | R²Y and Q² are close; record the gap | R²Y of 0.9 with Q² of 0.2 is a memorised training set |
| `leakage` | scaling, filtering, feature selection, and sPLS-DA tuning ran inside each training fold; replicate measurements of one subject stay in one fold | Leakage inflates Q² and can pass a careless permutation test |
| `ncomp_selection` | the number of components was chosen by CV, and the rule is stated | Each extra component raises R²Y whether or not it predicts |

## Judgments and decisions

| Belongs to code | Belongs to the model |
| :--- | :--- |
| permutation p, Q², and gap against thresholds | Do the high-VIP features higher in one class share a pathway or chemical class? |
| which features pass VIP and stability cutoffs | Is the pattern consistent with the condition, or with a confounder (fasting state, medication, sample handling, haemolysis)? |
| direction of each feature | Which features are likely the same underlying signal (adducts, isotopes, fragments of one compound)? |

## Pitfalls

- Presenting the scores plot as evidence of separation.
- Treating VIP > 1 as a significance threshold. It marks above-average contribution; about a third
  of features exceed it in most models.
- Unbalanced classes: overall accuracy rewards predicting the majority. Use balanced error or AUC.
- In untargeted metabolomics, several features are one compound. Count compounds, not features,
  before calling a pathway enriched.
