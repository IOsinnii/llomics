---
name: llomics-limma
description: >
  Build an llomics payload from a feature-wise differential analysis (limma, limma-voom, and the
  same-shaped tables from DESeq2 and edgeR): one state per contrast with hit counts, top up and
  down features, a control panel, the shape of the p-value distribution, and trust checks. Use
  when a topTable or results table has thousands of rows and the question is what changed, how
  much, and whether the result can be believed.
---

# llomics limma

The full result is one row per feature per contrast, typically 15,000 to 60,000 rows. The decision
is per contrast: is there a real response, what is it, and which features carry it. Read
`llomics-core` first.

## Unit

One contrast. State its definition in words, including the direction:
`"treated − control, adjusted for batch; positive logFC = higher in treated"`. A reader who has to
guess the sign of the contrast will guess wrong half the time.

A second level of unit is useful when a contrast has hundreds of hits: one unit per enriched gene
set or per co-expression module among the hits, each with its member features.

## Extract

| Field | limma | DESeq2 | edgeR |
| :--- | :--- | :--- | :--- |
| full table | `topTable(fit, coef, number = Inf, sort.by = "none")` | `results(dds, contrast)` | `topTags(lrt, n = Inf)$table` |
| effect | `logFC` | `log2FoldChange` (say if shrunk, and how) | `logFC` |
| abundance | `AveExpr` | `baseMean` | `logCPM` |
| stat | `t` | `stat` | `LR` or `F` |
| adjusted p | `adj.P.Val` | `padj` | `FDR` |
| hit counts | `summary(decideTests(fit))` | count from `padj` | `summary(decideTests(lrt))` |
| variance prior | `fit$df.prior`, `fit$s2.prior` | dispersion trend | `common.dispersion`, prior df |

Write the full table to disk and record its path in `run.source_tables`.

## Reference extractor

`scripts/llomics_limma.R` implements this skill for a limma fit: `llomics_limma(fit, coefs, ...)`
writes the payload and the full table per contrast, and computes the trust checks below from the
fit and the sample table. `references/questions.contrast.json` is a starter set of questions and
gates for the `llomics-judge` ask step. DESeq2 and edgeR have no extractor yet; follow the table
above.

## Unit state

```json
{
  "unit_id": "treated_vs_control",
  "unit_type": "contrast",
  "definition": "treated − control, adjusted for batch; positive logFC = higher in treated",
  "size": {"n_treated": 6, "n_control": 6, "n_tested": 15234},
  "facts": {
    "fdr": 0.05,
    "n_up": 412,
    "n_down": 318,
    "n_up_lfc1": 96,
    "n_down_lfc1": 41,
    "n_near_miss": 207,
    "pi0": 0.74,
    "p_histogram": [3950, 1720, 1410, 1290, 1220, 1180, 1150, 1120, 1100, 1094],
    "logfc_quantiles": {"q01": -1.9, "q25": -0.21, "q50": 0.01, "q75": 0.23, "q99": 2.4}
  },
  "top_features": {
    "up": [{"feature": "ENSG00000185745", "symbol": "IFIT1", "effect": 4.21, "ave_expr": 6.3, "stat": 14.2, "p_adj": 3.1e-9}],
    "down": [{"feature": "ENSG00000170345", "symbol": "FOS", "effect": -2.87, "ave_expr": 7.9, "stat": -11.6, "p_adj": 4.4e-8}]
  },
  "panel": {
    "ISG15": {"effect": 3.9, "p_adj": 6.0e-9},
    "ACTB": {"effect": 0.02, "p_adj": 0.97},
    "XIST": null
  },
  "tail": {"n_shown": 30, "n_eligible": 730, "rule": "adj.P.Val < 0.05, ranked by |t|, 15 per direction"}
}
```

Compute in code and store in `facts`:

- **Counts at more than one threshold.** Hits at the FDR cutoff, hits that also pass an effect
  size cutoff, and `n_near_miss` (features between the cutoff and twice the cutoff). A contrast
  with 730 hits and 207 near misses reads differently from one with 730 hits and a cliff.
- **`p_histogram`**: raw p-values in ten equal bins. This one row tells more about the health of
  the analysis than any single gene.
- **`pi0`**: estimated share of true nulls, from the flat right side of the histogram.
- **`logfc_quantiles`**: where the bulk of effects sit, so a logFC of 1 can be placed.

Rank `top_features` by the moderated statistic, not by adjusted p-value: ties in adjusted p are
common and hide the order. Keep up and down separate so a strong direction does not crowd out the
other.

`panel` holds expected responders (known positives for the treatment), expected non-responders
(housekeeping), and technical flags (sex genes when sex is unbalanced). A feature removed by the
expression filter is `null`, with the filter listed in `run.filters`.

## Required trust checks

| Check | Pass when | Why it matters |
| :--- | :--- | :--- |
| `replicates` | each group has enough independent biological replicates; record the minimum n | With n = 2 the moderated variance is mostly prior. Technical replicates are not replicates |
| `pvalue_distribution` | the histogram is flat with a peak near 0 | A peak near 1, a hill in the middle, or a U shape means the model is wrong (missing covariate, bad variance model, correlated samples), and the hit list is not trustworthy at any cutoff |
| `batch_confounding` | batch is in the design and is not aliased with the contrast; record the group × batch table | If every treated sample is in one batch, the contrast is the batch |
| `multiple_testing` | method and scope are stated: per contrast or across all contrasts | Ten contrasts each at FDR 5% is not a 5% family |
| `mean_variance` | the mean–variance trend was modelled (voom, `trend = TRUE`, or a count model) | Otherwise low-abundance features dominate or vanish for technical reasons |

Also record when present: `duplicateCorrelation` for repeated measures, `treat` thresholds (the
p-values then test against a non-zero fold change), sample weights, and for DESeq2 the number of
features with `padj = NA` and why (independent filtering or Cook's distance outliers). An `NA` is
`null` in the payload, not 1.

## Judgments and decisions

| Belongs to code | Belongs to the model |
| :--- | :--- |
| hit counts, near misses, histogram shape tests | Does the up set form a recognisable program? Which one? |
| whether expected positives in the panel moved | Is the down set the inverse of the up set's program, or a separate one? |
| whether a contrast has too few hits to interpret | Are the top features consistent with the stated treatment, or with an artefact (stress, sex, cell-type shift)? |
| enrichment statistics (run a proper test) | Which enriched sets describe the same biology and should be read together? |

Do not ask a model to do enrichment by reading a gene list. Run the test in code, and ask the
model to interpret and group the result.

## Pitfalls

- A large logFC at low `ave_expr` is often a ratio of two near-zero numbers. Always send abundance
  next to effect.
- Sorting by adjusted p-value and taking the top 30 returns one direction when the response is
  asymmetric.
- "No significant genes" with a flat p-value histogram is a clean negative. With a peak near 1 it
  is a broken model. The payload must let the reader tell these apart.
