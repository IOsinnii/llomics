---
name: llomics-synthesis
license: MIT
description: >
  Combine llomics payloads from several methods run on the same data (for example limma, lasso,
  PLS-DA, and PCA) into one evidence table per feature, and separate agreement that strengthens a
  finding from agreement that only repeats it. Use when the same dataset was analysed more than
  one way and the question is which findings hold across methods and which belong to one method.
---

# llomics synthesis

Each method payload answers its own question. A feature that is a limma hit, is stably selected by
the lasso, and has a high VIP is better supported than one that appears once, but only if the
methods could have disagreed. Read `llomics-core` first.

## Before joining

1. **Same observations.** Compare `dataset.n_observations` and the sample lists. Methods run on
   different subsets do not corroborate each other.
2. **Same feature IDs.** Join on the stable identifier, never on symbols or names.
3. **Same direction convention.** Convert every method's direction to "higher in group X" before
   comparing. Loading and coefficient signs are not comparable across methods.
4. **Same preprocessing lineage.** If all methods share one normalised matrix, they share its
   artefacts. Record this; it bounds how independent the agreement is.

## Build the evidence table in code

Take the union of features that any payload kept in `top_features` or `panel`. For each, look up
its values in every method's full table through `run.source_tables`. Do not infer "not a hit" from
absence in a payload: the payload shows the top k, and a feature ranked k+1 is absent for a
different reason than one with no effect.

```json
{
  "feature": "ENSG00000185745",
  "symbol": "IFIT1",
  "higher_in": "treated",
  "evidence": {
    "limma": {"effect": 4.21, "p_adj": 3.1e-9, "rank": 1, "kept": true},
    "lasso": {"effect": 0.38, "selection_freq": 0.94, "kept": true},
    "plsda": {"vip": 2.3, "higher_in": "treated", "kept": true},
    "pca": {"component": "PC1", "loading_rank": 4, "kept": true}
  },
  "n_methods_supporting": 4,
  "direction_consistent": true,
  "tested_in": ["limma", "lasso", "plsda", "pca"]
}
```

`tested_in` matters: a feature filtered out before the lasso was not rejected by it. Use `null`
for a method that never saw the feature.

The unit of the synthesis payload is the feature, or better, a group of features that the methods
agree on. Its run card lists the source payloads with their hashes, and its trust layer carries
forward every `warn`, `fail`, and `not_run` from them.

## What agreement means

| Pattern | Reading |
| :--- | :--- |
| Univariate hit and multivariate selection, same direction | Strongest support: the feature differs on its own and adds to prediction |
| Univariate hit, never selected by lasso, has a selected stand-in | Same signal, represented by a correlated feature. Count the group once |
| High VIP or selected, no univariate difference | Contributes through combination with other features. Plausible, harder to validate singly |
| Hit in limma, top loading on a component associated with batch | Suspect: the contrast may be partly the batch |
| Direction differs between methods | Check the direction conventions first. If real, suspect a confounder or an interaction |
| Agreement among limma, lasso, and PLS-DA on the same samples and outcome | Expected to a large degree: all three see the same group difference. It confirms robustness to method, not replication |

Only an independent cohort, a held-out set, or an orthogonal assay is replication. Say so in the
payload when none exists.

## Judgments and decisions

| Belongs to code | Belongs to the model |
| :--- | :--- |
| the join, the counts, direction consistency | Do the features supported across methods form one program or several? |
| grouping features by correlation | Which single-method findings are plausible biology that one method is better placed to see? |
| carrying trust statuses forward | Which conflicts between methods point to a confounder, and which one? |

## Pitfalls

- Counting methods as votes. Three methods that share one normalised matrix and one outcome are
  closer to one vote than three.
- Intersecting top-k lists. The intersection shrinks with k and with noise in the ranking; use
  the full tables.
- Dropping the trust layer during the join. A synthesis built on a PLS-DA with no permutation
  test inherits that gap.
