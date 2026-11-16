# llomics

Skills that teach a coding agent to turn the output of statistical and ML steps in omics pipelines
into small, self-describing payloads that an LLM can reason over.

A single script run can produce a 20,000-row differential expression table, a loadings matrix, or
a regularisation path. No person reads those whole, and pasting them into a prompt buries the
signal. llomics defines what to send instead: one screen per unit of judgment, with the context and
the caveats attached.

Status: v0.1, a first draft of the contract and nine skills. One worked example on real data
(`examples/pbmc10k`). limma has a tested R extractor and a simulated example
(`examples/limma_sim`). The skills for PCA, GLM, lasso, and PLS-DA have not yet been exercised on
a real pipeline, and the Bedrock call in the ask step has not been run against a live endpoint.

## The pattern

```text
pipeline object  ──extract──▶  payload (JSON)  ──ask──▶  typed answers  ──gate──▶  auto / review
 415 MB AnnData                 19 units, ~1.6 KB each     per unit              reason codes
 (L2, stays on disk)            (L0 run card + L1 units)   (model)               (code)
```

| Layer | What it is | Who reads it |
| :--- | :--- | :--- |
| L0 run card | what was run, the filters, the checks that passed or were not run | model and person |
| L1 unit states | one small record per cluster, component, contrast, or model | model, one unit per request |
| L2 source tables | the full result, referenced by path | code, on demand |

## Skills

| Skill | Covers |
| :--- | :--- |
| [`llomics-core`](skills/llomics-core/SKILL.md) | The payload contract, the nine rules, the trust layer, the validator |
| [`llomics-clusters`](skills/llomics-clusters/SKILL.md) | Clustering plus marker ranking (Scanpy, Seurat) |
| [`llomics-pca`](skills/llomics-pca/SKILL.md) | PCA components, loadings, metadata association |
| [`llomics-limma`](skills/llomics-limma/SKILL.md) | Feature-wise differential analysis (limma, DESeq2, edgeR) |
| [`llomics-glm`](skills/llomics-glm/SKILL.md) | A fitted GLM as the unit: coefficients and diagnostics |
| [`llomics-lasso`](skills/llomics-lasso/SKILL.md) | Lasso and elastic net: selection, stability, stand-ins |
| [`llomics-plsda`](skills/llomics-plsda/SKILL.md) | PLS-DA, sPLS-DA, OPLS-DA: VIP, cross-validation, permutation |
| [`llomics-judge`](skills/llomics-judge/SKILL.md) | Typed questions per unit, gates in code, reason codes |
| [`llomics-synthesis`](skills/llomics-synthesis/SKILL.md) | Joining payloads from several methods on the same data |

Every method skill has the same sections: the unit, where each field comes from in R and Python,
an example unit state, the required trust checks, what belongs to code and what to the model, and
the common misreadings.

## Install

The repository follows the `skills/<name>/SKILL.md` layout used by the `skills` CLI:

```bash
npx skills add <owner>/llomics
```

Or copy the folders under `skills/` into your agent's skills directory.

## Validate a payload

```bash
python skills/llomics-core/scripts/validate_payload.py examples/pbmc10k/payload.json
```

Standard library only. Errors (missing trust checks, zero p-values, wrong tail counts, NaN) exit
non-zero.

## Adapting it to your stack

The payload is plain JSON, and it sits between the analysis and the model. Neither side knows
about the other, so each can be replaced on its own:

```text
analysis framework ──▶ extractor ──▶ payload.json ──▶ ask step ──▶ gates in code
(swap freely)          (rewrite)     (unchanged)      (rewrite)     (unchanged)
```

### A different analysis framework (Seurat, DESeq2, mixOmics, your own code)

Only the extractor changes. The contract, the validator, the questions, and the gates stay as
they are.

1. Open the method skill and find its **Extract** table. It maps each payload field to the
   object or column that holds it in the common R and Python tools (for example `effect` is
   `logfoldchanges` in Scanpy and `avg_log2FC` in Seurat).
2. Write the extractor in the language your pipeline already uses. It only has to emit JSON.
3. Run the validator on the result. It reads the JSON file and does not care what produced it.

The same step from Seurat, in R:

```r
markers <- FindAllMarkers(obj, only.pos = TRUE)
top_features <- lapply(split(markers, markers$cluster), function(m) {
  m <- head(m[order(m$p_val_adj, -m$avg_log2FC), ], 15)
  data.frame(feature = m$gene, effect = m$avg_log2FC,
             p_adj = pmax(m$p_val_adj, 1e-300),   # clamp underflow to run.p_floor
             pct_in = m$pct.1, pct_out = m$pct.2)
})
# ...assemble run, dataset, trust, and units as in llomics-core, then:
jsonlite::write_json(payload, "payload.json", auto_unbox = TRUE, null = "null", na = "null", digits = NA)
```

If your tool is not in the table, add a column: decide which of its outputs are the `effect`, the
`stat`, and the `p_adj`, and record the tool and version in `run.tool`.

For a method llomics does not cover yet (WGCNA modules, random forest importance, a survival
model), write a new skill with the same sections: choose the unit of judgment, define its state,
list the trust checks, and add the check names under a new key in `REQUIRED_CHECKS` in
`skills/llomics-core/scripts/validate_payload.py`.

### A different model (Claude on Amazon Bedrock, the Anthropic API, any other LLM)

Only the ask step changes. The payload, the truth table, and the gates stay as they are. The
worked example uses TypeSafe Jev; any model that can return answers in a fixed schema works. Map
the three question shapes from `llomics-judge` onto a schema:

| Question shape | In a JSON schema |
| :--- | :--- |
| Choice among named options | a string restricted to the option names |
| Is this statement true? | a boolean |
| Score on described levels | a string restricted to the level descriptions |

With Claude on Amazon Bedrock, through the Anthropic Python SDK and structured outputs:

```python
import json
from typing import Literal

from anthropic import AnthropicBedrockMantle
from pydantic import BaseModel


class ClusterJudgment(BaseModel):
    cell_type: Literal["CD4+ T cells", "CD8+ T cells", "NK cells", "B cells", "CD14+ Monocytes"]
    any_listed_type_fits: bool
    markers_define_one_lineage: bool
    exclusive_lineages_coexpressed: bool


client = AnthropicBedrockMantle(aws_region="us-east-1")  # uses your AWS credentials


def judge(unit: dict, reference: dict) -> ClusterJudgment | None:
    response = client.messages.parse(
        model="anthropic.claude-opus-5-5",
        max_tokens=16000,
        system="You judge one cluster from a single-cell RNA-seq analysis. Use only the state provided.",
        messages=[{"role": "user", "content": json.dumps({"unit": unit, "dataset_reference": reference})}],
        output_format=ClusterJudgment,
    )
    if response.stop_reason == "refusal":
        return None  # no answer: send the unit to manual review, never treat it as a "no"
    return response.parsed_output
```

A complete version of this step, with repeated runs, vote shares, and gates, is
`skills/llomics-judge/scripts/judge_bedrock.py` (see the HPC section below).

For the Anthropic API directly, the client is `anthropic.Anthropic()` and the model ID is
`claude-opus-5-5`; the rest is identical. For another provider, replace the client call and keep
the schema.

Three things change when you leave a typed-judgment model, and the gates must account for them:

- **No calibrated probabilities.** Jev returns a probability with each answer. A general LLM
  returns the answer only, and a confidence number it reports about itself is not calibrated. To
  get a usable number, send the same request several times and use the share of runs that agree
  (5 runs with 4 agreeing gives 0.8). Gate on that share and on the margin between the top two
  answers.
- **Thresholds do not transfer.** A cutoff of 0.85 set for one model's probabilities means
  something else for another model's vote share. Set the thresholds again, and write the model ID
  next to them in the output.
- **Trust gates and fact gates are unaffected.** They read the payload, not the model, which is
  one more reason to keep measured facts out of the questions.

## Running it on an HPC with Amazon Bedrock

For a cluster where the only model access is Claude through Amazon Bedrock. The analysis and the
extractor are R; everything after the payload is Python.

**1. Requirements.** R with `limma` and `jsonlite`. Python 3.10 or later. The validator and the
dry run need no packages; the real ask step needs `pip install anthropic`.

**2. Smoke test, no model involved.** From the repository root:

```bash
Rscript examples/limma_sim/run_example.R
```

```bash
python skills/llomics-core/scripts/validate_payload.py examples/limma_sim/out/limma-sim.llomics.json
```

```bash
python skills/llomics-judge/scripts/judge_bedrock.py examples/limma_sim/out/limma-sim.llomics.json skills/llomics-limma/references/questions.contrast.json --dry-run
```

The dry run prints the exact request that would be sent for the first unit and evaluates the
gates that do not need the model. Nothing leaves the machine.

**3. Your own fit.** After `eBayes()` in your pipeline:

```r
source("skills/llomics-limma/scripts/llomics_limma.R")
llomics_limma(
  fit, coefs = "treated_vs_control",
  definitions = c(treated_vs_control = "treated - control, adjusted for batch; positive logFC = higher in treated"),
  targets = targets, group = "group", batch = "batch", batch_in_design = TRUE,
  analysis_id = "my-study", question = "Which genes respond to treatment?",
  steps = c("filterByExpr", "TMM", "voom", "lmFit ~0 + group + batch", "eBayes"),
  mean_variance = "voom", panel = c("ISG15", "ACTB"), symbol_col = "symbol",
  out_dir = "outputs"
)
```

It writes `outputs/my-study.llomics.json` and the full table per contrast, and computes the trust
checks it can from the fit and the sample table. What it cannot know (whether batch is in the
design, whether voom was used) you declare; anything left undeclared is recorded as `not_run`.

**4. Ask Claude on Bedrock.** One unit, one run first:

```bash
python skills/llomics-judge/scripts/judge_bedrock.py outputs/my-study.llomics.json skills/llomics-limma/references/questions.contrast.json --region us-east-1 --runs 1 --limit 1
```

Then the full run (three requests per unit by default, so that agreement between runs can be
measured):

```bash
python skills/llomics-judge/scripts/judge_bedrock.py outputs/my-study.llomics.json skills/llomics-limma/references/questions.contrast.json --region us-east-1 --out outputs/my-study.answers.json
```

`answers.json` records the payload hash, the model, every run, the vote shares, and for each unit
`auto` or `review` with reason codes.

Things that differ between clusters:

- **Model ID.** The default is `anthropic.claude-opus-5-5`. Use `--model` or `LLOMICS_MODEL` for
  whichever Claude model your account has enabled in that region.
- **Endpoint.** The default client uses Bedrock's Messages-API endpoint. If your organisation only
  allows the older `bedrock-runtime` InvokeModel path, add `--client legacy`.
- **Network.** Compute nodes often have no outbound internet. Run the extractor wherever the
  analysis runs, and the ask step wherever HTTPS to Bedrock is allowed (a login node, or through
  your site's `HTTPS_PROXY`). It is a handful of small requests and needs no job allocation.
- **Credentials.** The SDK reads the standard AWS chain: environment variables, a shared profile,
  or an attached role. Do not put keys in scripts or job files.
- **What leaves the cluster.** Only the request the dry run prints: the run description, hit
  counts, group names and sizes, and statistics for the top and panel features. No per-sample
  values. Confirm this is within your data policy before the first real request.

`questions.contrast.json` holds the questions and the gate thresholds. It is a starter written for
human or mouse transcriptomics; edit it for your data.

## Worked example: PBMC 10k

`examples/pbmc10k` converts the cluster summaries of the
`pbmc10k-jev-annotation` pipeline into the contract. That pipeline is where
the pattern comes from: Scanpy produces a 415 MB AnnData, a summary step reduces it to 19 cluster
states, TypeSafe Jev answers four typed questions per cluster, and five gates in code decide
between auto-assignment and manual review.

```bash
python examples/pbmc10k/to_llomics.py path/to/cluster_summaries.json examples/pbmc10k/payload.json
```

What the case study already did, which the contract keeps:

- One state per cluster, small enough to judge at a glance.
- Dataset-wide reference values sent with each cluster, so QC is relative.
- A data-driven marker list next to a fixed canonical panel.
- Unmeasured genes recorded as `null`, not `0.0`.
- Marker lists never padded with non-significant genes.
- Judgments from the model, decisions in code, with a reason code for every failed gate.

What converting it showed, which the contract adds:

| Finding in the case study | Contract rule |
| :--- | :--- |
| 71 of 285 marker rows carry `pval_adj: 0.0` | Clamp to a declared floor; rank by the statistic |
| The pipeline parameters and the marker exclusion pattern live only in the script | Run card with `steps` and `filters` |
| No record of how many markers passed the filter before the top 15 were taken | `tail.n_eligible` |
| One donor, cells as replicates, p-values near 1e-300 | Required `pseudoreplication` check |
| The doublet question's answers followed the Scrublet fraction; `pct_predicted_doublet >= 5` reproduces the gate on all 19 clusters | Measured facts are gated in code; the model is asked only for the judgment |
| The reference is a median of cluster medians, computed in the annotation step | Reference computed at extraction, stored in the payload with its definition |
| Marker rows have fold change but no fraction of cells expressing | `pct_in`, `pct_out` in `llomics-clusters` |

## Open questions

- **Schema strictness.** The validator enforces the envelope and the trust checks, not the
  method-specific `facts`. A JSON Schema per method would catch more and constrain more.
- **Feature-level units.** When the decision really is per feature (annotating 700 hits one at a
  time), the unit state needs gene-level context from outside the pipeline. Not covered yet.
- **Reference extractors.** Only limma has one (`skills/llomics-limma/scripts/llomics_limma.R`).
  DESeq2, PCA, GLM, lasso, and PLS-DA are described in their skills but not yet implemented.
- **Refusals.** When Claude declines a request, the unit goes to review. Automatic fallback to a
  second model is not set up.
- **Thresholds.** Trust checks say "pass when" in words. Whether llomics should ship default
  numeric thresholds, or leave them to each project, is undecided.
