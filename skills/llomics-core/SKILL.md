---
name: llomics-core
description: >
  The llomics payload contract: how to turn the output of a statistical or ML step in an omics
  pipeline (PCA, limma/DESeq2, GLM, lasso, PLS-DA, clustering and marker ranking) into a compact,
  self-describing JSON payload that an LLM or a typed-judgment model can reason over. Use when a
  pipeline step produces more rows than a person can read, when writing or reviewing the
  "summary for the model" step of a pipeline, and before using any other llomics-* skill.
---

# llomics core: the payload contract

A pipeline step produces thousands of rows. A model can judge about one screen of them well. The
payload is that screen: the smallest state from which a domain expert would reach the same
conclusion as someone holding the full result.

The method skills (`llomics-pca`, `llomics-limma`, `llomics-glm`, `llomics-lasso`,
`llomics-plsda`, `llomics-clusters`) say what goes on the screen for each method. This skill
defines what they share.

## Three layers

| Layer | What it is | Size | Who reads it |
| :--- | :--- | :--- | :--- |
| L0 run card | What was run, on what, with which filters, and which checks passed | one screen | model and person, first |
| L1 unit states | One small record per unit of judgment | 1–2k tokens each | model, one unit per request |
| L2 source tables | The full result on disk, referenced by path | unlimited | code only, on demand |

The payload holds L0 and L1 and points to L2. Nothing in L2 is pasted into a prompt. When a
question needs a row that L1 dropped, code fetches it from the path in `run.source_tables`.

## Choose the unit of judgment

The unit is the thing a decision is made about. It is rarely the feature: one gene's row cannot be
judged alone, and there are too many.

| Method | Unit | Features appear as |
| :--- | :--- | :--- |
| Clustering + markers | cluster | markers of the cluster |
| PCA | component | loadings of the component |
| limma / DESeq2 / edgeR | contrast | hits of the contrast |
| GLM (one model) | model | terms of the model |
| lasso / elastic net | model at the chosen λ | selected features |
| PLS-DA | model, then component | high-VIP features |

Check each unit with the **human test**: could an expert, shown only this unit's state and the
dataset reference, answer the question at a glance? If not, the state is missing something or the
unit is wrong.

## The envelope

```json
{
  "llomics": "0.1",
  "method": "limma",
  "run": {
    "analysis_id": "stim-vs-ctrl-2026-11",
    "question": "Which genes respond to stimulation, and is the response an interferon program?",
    "tool": "limma 3.60 (voom, eBayes robust=TRUE)",
    "steps": ["filterByExpr", "TMM", "voom", "lmFit ~0 + group + batch", "eBayes"],
    "filters": [{"rule": "adj.P.Val < 0.05", "why": "FDR control"}],
    "p_floor": 1e-300,
    "source_tables": {"top_table": "outputs/limma_stim_vs_ctrl.tsv"}
  },
  "dataset": {
    "n_observations": 12,
    "observation_type": "sample",
    "n_features": 15234,
    "reference": {}
  },
  "trust": [
    {"check": "replicates", "status": "pass", "value": 6, "note": "6 per group"}
  ],
  "units": [
    {
      "unit_id": "stim_vs_ctrl",
      "unit_type": "contrast",
      "facts": {},
      "top_features": [],
      "panel": {},
      "tail": {"n_shown": 30, "n_eligible": 730, "rule": "adj.P.Val < 0.05, ranked by |t|"}
    }
  ]
}
```

`top_features` is a list of rows, or an object of lists when direction matters
(`{"up": [...], "down": [...]}`). Every row has `feature` and uses the same column names across
methods where they apply: `effect`, `stat`, `p_adj`.

## Rules

1. **Two views of the features.** `top_features` is data-driven: the top k by the method's own
   ranking. `panel` is hypothesis-driven: a fixed list the analyst expects to matter (lineage
   markers, known positives, housekeeping and technical genes), reported whether or not they
   ranked. The top list shows what the data says. The panel shows what is absent, which a top list
   cannot.
2. **Null is not zero.** A panel gene that was filtered out or never measured is `null`, never
   `0.0`. A check that was not performed has status `not_run`, never a missing entry.
3. **The reference travels with the unit.** Put the dataset-wide values each unit is compared
   against in `dataset.reference`, and send them with every unit. Then compute the comparison in
   code and store the result in `facts` (`"doublet_fraction_vs_reference": 5.5`). The model should
   read a comparison, not perform one.
4. **Never pad, always count the tail.** If only 9 features pass the filter, send 9. `tail`
   records how many were eligible and the rule that ranked them, so the reader knows whether 15
   rows are all of the signal or 2% of it.
5. **Declare every filter.** Any rule that removed features before ranking (expression filter,
   excluded gene patterns, significance and direction cutoffs) goes in `run.filters` with its
   reason. An undeclared filter looks like a biological absence.
6. **No zero p-values, and say what p means.** A p-value of exactly 0 is an underflow: clamp it to
   `run.p_floor` and keep the test statistic, which still ranks. If the test's replicates are not
   independent (cells from one donor, technical replicates), record that in `trust` so the
   p-values are read as a ranking and not as evidence.
7. **Stable ID plus readable name.** Keep the identifier that joins to annotation
   (`ENSG…`, UniProt, HMDB) and add the symbol or name a person would recognise.
8. **One unit, one screen.** Aim for 1–2k tokens per unit. If a unit needs more, the unit is too
   coarse, or columns are being sent that no question reads.
9. **The payload is a file.** Write it to disk next to the other processed intermediates, commit
   it when it is small, and record its hash with any model answers derived from it. Model answers
   vary between runs; the payload must not.

## The trust layer

`trust` is a list of named checks. Each method skill lists its required checks. Every required
check must be present with one of five statuses:

| Status | Meaning |
| :--- | :--- |
| `pass` | Performed, no concern |
| `warn` | Performed, concern noted in `note` |
| `fail` | Performed, result should not be interpreted as it stands |
| `not_run` | Not performed. The reader must assume the risk is live |
| `not_applicable` | Does not apply to this design; say why in `note` |

The extractor never refuses to write a payload. It states what is missing, and the code that
consumes the answers decides what a `fail` or `not_run` blocks (see `llomics-judge`).

## Validate

```bash
python skills/llomics-core/scripts/validate_payload.py path/to/payload.json
```

The validator (standard library only) checks the envelope, the required trust checks for the
method, zero or out-of-range p-values, NaN and infinity, tail counts, and the per-unit size
budget. Errors exit non-zero. Warnings list `not_run` checks and unknown tail sizes.

## After the payload

- `llomics-judge`: ask typed questions of each unit and turn answers into decisions in code.
- `llomics-synthesis`: join payloads from several methods run on the same data.
