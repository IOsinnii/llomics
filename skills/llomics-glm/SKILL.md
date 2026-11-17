---
name: llomics-glm
license: MIT
description: >
  Build an llomics payload from a generalised linear model (R glm, statsmodels GLM or Logit,
  logistic, Poisson, negative binomial): one state per fitted model with coefficients on an
  interpretable scale, fit statistics, and diagnostics for convergence, separation,
  overdispersion, and collinearity. Use when a model has many terms or many models were fitted,
  and the question is which effects are real and whether the fit can be trusted.
---

# llomics GLM

A single GLM summary is short, but the things that decide whether it can be believed (dispersion,
separation, leverage, collinearity) are not in the coefficient table. With many terms, or one model
per outcome, the tables stop being readable. Read `llomics-core` first.

For one GLM per feature (DESeq2, edgeR, per-gene logistic models) the unit is the contrast and the
state follows `llomics-limma`. This skill covers the case where the model is the unit.

## Unit

One fitted model. Its "features" are its terms.

## Extract

| Field | R | statsmodels |
| :--- | :--- | :--- |
| coefficients | `summary(fit)$coefficients` | `res.params`, `res.bse`, `res.tvalues`, `res.pvalues` |
| intervals | `confint(fit)` (profile) or `confint.default(fit)` (Wald) | `res.conf_int()` |
| deviance | `deviance(fit)`, `fit$null.deviance` | `res.deviance`, `res.null_deviance` |
| information criteria | `AIC(fit)` | `res.aic` |
| dispersion | `sum(residuals(fit, "pearson")^2) / df.residual(fit)` | `res.pearson_chi2 / res.df_resid` |
| convergence | `fit$converged`, `fit$iter` | `res.converged` |
| collinearity | `car::vif(fit)` | `variance_inflation_factor` |
| influence | `cooks.distance(fit)` | `res.get_influence().cooks_distance[0]` |

## Unit state

```json
{
  "unit_id": "response_logit",
  "unit_type": "model",
  "definition": "logit P(responder) ~ age + sex + stage + log2(CRP) + treatment; reference levels: sex=F, stage=I, treatment=placebo",
  "size": {"n": 412, "n_events": 97, "n_dropped_missing": 23, "n_terms": 7},
  "facts": {
    "family": "binomial",
    "link": "logit",
    "effect_scale": "odds ratio per unit, or versus reference level",
    "null_deviance": 449.6,
    "residual_deviance": 391.2,
    "df_residual": 405,
    "aic": 405.2,
    "mcfadden_r2": 0.13,
    "dispersion": 1.04,
    "events_per_term": 13.9,
    "max_vif": 2.1,
    "n_high_influence": 3
  },
  "top_features": [
    {"feature": "treatment[active]", "effect": 2.84, "ci_low": 1.71, "ci_high": 4.72, "coef": 1.044, "se": 0.259, "stat": 4.03, "p": 5.6e-5},
    {"feature": "log2_crp", "effect": 0.78, "ci_low": 0.66, "ci_high": 0.92, "coef": -0.248, "se": 0.085, "stat": -2.92, "p": 0.0035}
  ],
  "panel": {},
  "tail": {"n_shown": 2, "n_eligible": 7, "rule": "all terms when 30 or fewer; otherwise top by |z| plus prespecified terms in panel"}
}
```

- **Report the effect on the scale people reason in.** `effect` is the exponentiated coefficient
  for log and logit links (rate ratio, odds ratio), with its interval. Keep the raw `coef` and
  `se` too. State the scale once in `facts.effect_scale`.
- **Name the reference level** for every categorical term in `definition`. `sex[M]: 1.4` means
  nothing without it.
- **State units for continuous terms.** An odds ratio per year of age and per decade differ
  tenfold in the exponent.
- **Send all terms when there are 30 or fewer.** Above that, send the top terms by absolute z and
  put prespecified terms of interest in `panel`.
- **Flag interactions.** When a term appears in an interaction, its main effect is the effect at
  the reference level of the other term. Add `"conditional_on": "treatment=placebo"` to the row.

## Required trust checks

| Check | Pass when | Why it matters |
| :--- | :--- | :--- |
| `convergence` | the fitter converged within its iteration limit | Estimates from a non-converged fit are wherever the optimiser stopped |
| `separation` | no coefficient has an absurd magnitude with a huge standard error; no fitted probabilities at 0 or 1 | Under separation the estimate is infinite and the Wald p-value is meaningless. Use Firth or a penalised fit |
| `overdispersion` | Pearson dispersion is near 1 for Poisson and binomial counts | Above about 1.5, standard errors are too small and p-values too optimistic. Use quasi-likelihood or negative binomial |
| `collinearity` | maximum VIF is modest (below 5 is a common line) | Collinear terms trade effect between them; individual coefficients are unstable though the fit is fine |
| `events_per_variable` | events per term is at least about 10 for binomial models (`not_applicable` otherwise) | Below that, estimates are biased away from zero and intervals are unreliable |

Also record: how missing data were handled and how many rows were dropped, whether terms were
chosen by a selection procedure (stepwise selection invalidates the p-values), and the
multiple-testing scope when the same model is fitted to many outcomes.

## Judgments and decisions

| Belongs to code | Belongs to the model |
| :--- | :--- |
| which intervals exclude the null | Is the direction of each effect consistent with prior knowledge? |
| every diagnostic threshold | Which significant terms are likely proxies for an unmeasured confounder? |
| ranking of terms | Does the set of effects tell one coherent story, and what would contradict it? |

## Pitfalls

- Comparing coefficient magnitudes across terms on different scales. Standardise continuous
  predictors or do not rank by magnitude.
- Reporting pseudo-R² as variance explained. McFadden values of 0.2 to 0.4 are a good fit.
- Reading a non-significant term as "no effect" without its interval. The interval says whether
  the data exclude a meaningful effect or cannot tell.
