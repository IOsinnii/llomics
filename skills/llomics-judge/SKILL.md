---
name: llomics-judge
license: MIT
description: >
  How to consume an llomics payload: split each decision into measured facts (code) and judgments
  (model), ask typed questions of one unit at a time, and turn the answers into auto-accept or
  manual-review decisions with reason codes. Use when writing the annotation, triage, or
  interpretation step that follows a payload, with a typed-judgment model such as TypeSafe Jev, a
  structured-output LLM, or a reasoning agent.
---

# llomics judge: from payload to decision

The payload (`llomics-core`) is the state. This skill covers what is asked of it and what code does
with the answers. The pattern follows the PBMC case study in `examples/pbmc10k`: 19 cluster
states, four typed questions each, five gates in code, 8 clusters auto-assigned and 11 sent to
review with reasons.

## Start from the truth table

Write down the decision first: what gets auto-accepted, what goes to a person, and what is
discarded. List the inputs that decision needs. Then sort each input:

| Input | Belongs to | Example |
| :--- | :--- | :--- |
| A number already in the payload | code | `pct_predicted_doublet >= 5` |
| A comparison to the reference | code, stored in `facts` | `genes_vs_reference < 0.5` |
| A trust check status | code | `permutation_test == "not_run"` |
| Recognition that needs domain knowledge | model | "Do these markers define one lineage?" |

Only the last row becomes a question.

**Why the split matters.** In the case study, the `is_doublet` question joined two conditions with
"and/or": co-expression of exclusive lineage markers, and a high Scrublet fraction. The answers
followed the number. A code rule `pct_predicted_doublet >= 5` reproduces the doublet gate on all
19 clusters, including cluster 14, which was scored 0.90 for doublet while also scoring the
highest lineage specificity in the run (1.88 of 2). The measured half belonged in code. The
judgment half (are mutually exclusive lineages co-expressed?) was never asked on its own.

## Ask one unit at a time

Send one unit per request: the unit's state, `dataset.reference`, and any run-level trust notes
the question depends on. Do not send the whole payload and ask about unit 7. Requests for
different units are independent, so run them concurrently.

Each question asks for one judgment a person could make at a glance. If a question contains "and"
or "or", split it. If the model is TypeSafe Jev, the `typesafe-ai` skill covers question wording,
criteria, and how to read probabilities; follow it.

Useful question shapes per unit:

| Shape | Use for | Example |
| :--- | :--- | :--- |
| Choice among named options | identity | cell type, pathway family, technical vs biological |
| Probability a statement is true | presence or absence | "the up-regulated set is an interferon response" |
| Score on described levels | degree | how cleanly the features define one program |

A Choice always picks something. Pair it with a separate question that asks whether anything fits.

## Gate in code

```python
def triage(unit, answers, trust, thresholds):
    reasons = []
    # 1. Trust gates come first and do not depend on the model.
    for check in BLOCKING_CHECKS:
        if trust[check] in {"fail", "not_run"}:
            reasons.append(f"trust:{check}:{trust[check]}")
    # 2. Fact gates read numbers from the payload.
    if unit["facts"]["doublet_fraction_vs_reference"] >= thresholds["doublet_ratio"]:
        reasons.append("doublet_fraction")
    # 3. Judgment gates read model answers.
    if answers["identity"].margin < thresholds["margin"]:
        reasons.append("small_margin")
    return ("review" if reasons else "auto"), reasons
```

- **Trust gates first.** Decide per method which checks are blocking. A PLS-DA model with no
  permutation test should not reach auto-interpretation however confident the answer is.
- **Every failed gate leaves a reason code.** The reviewer's queue is sorted by reason.
- **Thresholds are declared design choices.** Keep them in one place, write them into the output,
  and do not tune them on the units they gate.
- **Use the margin, not only the top probability.** Top-1 minus top-2 separates "confidently X"
  from "X or Y".

## Reference ask step for Claude on Amazon Bedrock

`scripts/judge_bedrock.py payload.json questions.json` sends one request per unit per run to
Claude on Bedrock with a JSON schema built from the questions, then applies the gates in the order
above. `--dry-run` prints the exact request and evaluates the trust and fact gates without calling
a model.

Claude does not return a calibrated probability with each answer, so the script asks each unit
several times (`--runs`, default 3) and gates on the share of runs that agree. A request that
fails, is refused, or is cut off is recorded as no answer and sends the unit to review. It is
never counted as a "no".

The questions file holds `questions` (each with `id`, `shape` of `choice`, `bool`, `score`, or
`text`, and `text`) and `gates`:

- `blocking_checks`: an object naming, for each trust check, the statuses that send the unit to
  review. A warning should block for some checks and not for others: in the starter file it blocks
  for `replicates` and `pvalue_distribution`, and not for `batch_confounding`.
- `facts`: comparisons on numbers in the unit's `facts`, such as `n_hits < 10` or
  `n_technical_flags > 0`. A fact that is `null` sends the unit to review as missing.
- `judgments`: for a question, the accepted answers and the minimum share of runs that must give
  one of them.

`text` answers are recorded for the reviewer and are not gated.

The gate settings were changed once already by the benchmark in the llomics repository: under a
blanket "block on fail or not_run" rule, an analysis with two replicates per group and one with
batch left out of the design both passed. Run `python benchmark/score.py` after changing them.

## Record enough to rerun

For each unit store: the payload hash, the exact request, the model identifier, every answer with
its probabilities, the gate results, and the thresholds. Model answers vary between runs. For
units within a small distance of a threshold, repeat the request a few times and report the
spread, or send them to review.

## Reasoning mode

When the consumer is a reasoning model writing an interpretation instead of typed answers:

1. Give it the run card and trust layer first, then the unit states.
2. Require each claim to cite a `unit_id` and the `feature` rows it rests on.
3. Require it to state which claims the `warn`, `fail`, and `not_run` checks weaken.
4. When it needs a row that is not in a unit state, it reads `run.source_tables` with code. It
   does not guess the value.
