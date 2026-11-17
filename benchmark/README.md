# Benchmark

Eleven simulated two-group limma analyses, each with a known truth and a stated correct decision:
let the interpretation through (`auto`) or send it to a person (`review`). The benchmark measures
whether the gates and the model's judgments catch what they should, and whether they leave sound
analyses alone.

The scenarios use real human gene symbols so that a model can judge the hit lists. Expression
values are simulated.

| File | Purpose |
| :--- | :--- |
| `scenarios.json` | The truth of each scenario, the expected decision, and which layer should catch it |
| `make_scenarios.R` | Regenerates the payloads (needs limma and the llomics R package) |
| `payloads/` | One payload per scenario, committed so the benchmark runs without R |
| `score.py` | Scores the gates alone, or the gates plus model answers |

## Run it

Gates that need no model:

```bash
python benchmark/score.py
```

With Claude on Amazon Bedrock (three requests per scenario by default; answers are saved under
`benchmark/answers/<model>/` and can be re-scored with `--answers`):

```bash
python benchmark/score.py --ask --region us-east-1
```

## Results so far

Gates alone, with the starter gate settings in
`skills/llomics-limma/references/questions.contrast.json`:

| Scenario | Expected | Should be caught by | Result without a model |
| :--- | :--- | :--- | :--- |
| `clean_interferon` | auto | | passes the gates |
| `clean_hypoxia` | auto | | passes the gates |
| `null` | review | fact | caught: `too_few_hits` |
| `batch_aliased` | review | trust | caught: `batch_confounding` fail |
| `no_batch_info` | review | trust | caught: `batch_confounding` not_run |
| `unmodelled_batch` | review | trust | caught: `pvalue_distribution` warn |
| `sex_imbalance` | review | fact | caught: `technical_signature` (sex panel) |
| `stress_handling` | review | fact | caught: `technical_signature` (stress panel) |
| `composition` | review | model | pending |
| `two_replicates` | review | trust | caught: `replicates` warn |
| `technical_replicates` | review | none | not caught (known blind spot) |

The model has not been run on the benchmark yet. Until it has, nothing is known about the
`composition` scenario, or about whether the model raises false alarms on the two clean ones.
`clean_hypoxia` is the one to watch: the starter questions list hypoxia genes under handling
stress, so the model has to anchor on the stated comparison.

## What the first run changed

Under the original gate settings, a check blocked only on `fail` or `not_run`. Two scenarios got
through: two replicates per group (`replicates` was `warn`) and a design with batch left out
(`pvalue_distribution` was `warn`). The gate settings now name the blocking statuses per check, and
a warning blocks for those two. `tests/test_benchmark.py` pins this.

The sex-imbalance and handling-stress scenarios were first left to the model. They are now caught
in code: the packages score six technical panels against every hit list, and the gate reads the
resulting flag. In both scenarios the panel holds over 85% of the top hits; in the nine others no
panel is flagged.

## Known blind spot

`technical_replicates`: six measurements of three biological samples per group, analysed as six
independent samples. The test is anticonservative and returns over a hundred false hits. Nothing
in a fit object or a sample table reveals this; the analyst has to declare which column identifies
the biological sample. A `replicate_of` argument to the extractors would close it.

## Limits

- One simulated dataset per scenario. The results say whether a mechanism works, not how often.
- Only limma contrasts. Cluster payloads have no benchmark yet.
- The planted artefacts are clean-cut. Real ones are mixed with real signal.
