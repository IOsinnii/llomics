# How llomics is built: a beginner's walkthrough

This file is for someone who can write a script but has not built a package before, and wants to
see how one comes together on a real example. Each section is a question I had as a beginner,
answered by pointing at files in this repository. The last section lists what still stands between
this repository and a release on PyPI or CRAN.

Nothing here is needed to *use* llomics. For that, read [README.md](README.md).

## What does the package do, in one paragraph?

A differential expression test returns a table with thousands of rows. llomics turns that table
into a small JSON file, the *payload*, that holds what a reader needs to judge the result: the top
features, a few expected features, hit counts, and quality checks. A language model can then be
asked questions about each part, and plain code decides what its answers are allowed to settle.

## Why did it start with a document and not with code?

The first commit has no package in it. It has `skills/llomics-core/SKILL.md`, which describes what
a payload must contain, and a validator script that checks a file against that description.

That order matters. Before writing any code that *produces* a payload, there was a written
definition of a correct one and a program that could say yes or no. Every later piece of code had
a target to hit.

The definition now lives in one data file, `contract/llomics-contract.json`. A piece of it:

```json
"dge": {
  "unit_types": ["contrast"],
  "required_checks": ["replicates", "pvalue_distribution", "batch_confounding",
                      "multiple_testing", "mean_variance"]
}
```

This says: a differential-analysis payload is made of contrasts, and it must report on these five
quality checks. The validator reads this file. So do both packages. Changing a rule means editing
one line of JSON, not hunting through code.

## Why are there two packages, and why is the R one not just a copy?

Omics pipelines are written in R and in Python, so there is a package for each: `python/llomics`
and `r/llomics`. They are not translations of each other line by line. They share three things:

1. The contract file above.
2. The same behaviour: given the same table, both must produce the same payload.
3. A test that proves point 2 (see "How do I know R and Python agree?").

## Why is the code that talks to scanpy so short?

Each package is split in two layers.

**The core** takes a plain table with standard column names (`feature`, `effect`, `p`, `p_adj`)
and builds the payload. It knows nothing about scanpy, Seurat or limma. In Python it uses only the
standard library, so it runs anywhere and is easy to test. It is in `python/llomics/core.py` and
`r/llomics/R/core.R`.

**An adapter** knows one tool. It runs the tool's test, renames the tool's columns to the standard
names, collects what it can about the run, and hands everything to the core. The whole column
mapping for scanpy is one line in `python/llomics/adapters/scanpy.py`:

```python
COLUMNS = {"names": "feature", "logfoldchanges": "effect", "scores": "stat", "pvals": "p",
           "pvals_adj": "p_adj", "pct_nz_group": "pct_in", "pct_nz_reference": "pct_out"}
```

The payoff: every rule about what goes into a payload exists once per language, in the core.
Supporting a new tool means writing a short adapter, not a new implementation. The limma adapter
is about 100 lines; the core is about 600.

The first version did not have this split. The scanpy wrapper did everything itself (commit
`9440a9d`). The next commit (`8fd8fc8`) pulled the shared logic out into the core. Writing it once
the simple way, seeing what repeats, then restructuring, is a normal order to work in.

## What happens when I call a wrapper?

Follow `llomics.rank_genes_groups(adata, "leiden", n_features=15, panel=[...])`:

1. The adapter calls `sc.tl.rank_genes_groups` for you.
2. It reads the result for each cluster and renames the columns.
3. It gathers per-cell metadata: QC values, doublet flags, batch, sample.
4. It calls `core.clusters_payload(...)`.
5. The core picks the eligible markers, ranks them, keeps the top 15, counts how many more were
   eligible, compares each cluster's QC with the whole dataset, runs the quality checks, and
   scores the technical gene panels.
6. The adapter writes the JSON file and returns the payload.

Steps 1 to 3 are tool-specific. Steps 4 and 5 are identical for every tool.

## How do I know it works?

By tests: small programs that run the code and compare the result with what it should be. Run them
all with:

```bash
python -m unittest discover -s tests
```

There are four kinds here, and each catches a different sort of mistake.

**Rule tests** check one behaviour each. From `tests/test_core.py`:

```python
def test_null_is_not_zero(self):
    self.assertIsNone(self.unit["panel"]["NOT_MEASURED"])
```

A gene that was not measured must come out as `null`, never as `0`. If someone later changes the
code and breaks that, this test fails and names the rule that broke.

**Golden-file tests** compare a whole payload with a saved, reviewed copy
(`tests/fixtures/expected_dge.json`). They catch changes nobody intended, anywhere in the output.
When a change *is* intended, the saved copy is regenerated with
`LLOMICS_UPDATE_GOLDEN=1 python -m unittest discover -s tests`, and the difference is read before
it is committed.

**Adapter tests** run the real tools on tiny data: a synthetic AnnData for scanpy, the
`pbmc_small` object for Seurat, a simulated fit for limma. They are skipped when the tool is not
installed, so the rest of the suite still runs on a bare machine.

**Parity tests** are described next.

## What is a fixture?

A fixture is input data made for testing. `tests/fixtures/make_fixtures.py` generates small tables
with awkward cases planted on purpose: a p-value that underflowed to zero, two rows tied on the
test statistic, a gene whose name should be excluded, a missing adjusted p-value. Real data might
contain none of these on the day you test. A fixture guarantees every one of them is exercised
every time.

The generator uses a fixed random seed, so it produces the same tables on every machine.

## How do I know R and Python agree?

`tests/test_parity.py` builds a payload from the fixture tables with the Python core, asks R to
build one from the same tables (`tests/parity/build_from_fixtures.R`), and compares the two field
by field. Text and counts must match exactly; decimals must match to about three digits, because
the two languages round a number like 1.2345 differently in the last place.

This is the test that lets two implementations exist without drifting apart. When technical
panels and enrichment were added (commit `270913c`), they were written in both languages and the
parity test confirmed zero differences before the commit.

## Why is the same file in four places?

`llomics-contract.json` exists in `contract/`, inside the Python package, inside the R package,
and inside a skill folder. That looks like a mistake, and normally it would be.

The reason is that each of those is installed on its own. Someone who installs only the R package
gets only the files inside `r/llomics`. So each needs its own copy.

The rule that keeps this safe: only the file in `contract/` is ever edited. `python tools/sync.py`
copies it to the other places, and a test (`test_copies_are_in_sync`) fails if any copy differs.
Duplication is acceptable when a program maintains it and a test guards it.

## How is the benchmark different from the tests?

Tests ask: does the code do what I wrote it to do? The benchmark asks a harder question: is what I
wrote it to do actually *right*?

`benchmark/` holds eleven simulated analyses where the truth is known because it was planted: one
has batch confounded with treatment, one compares all-female with all-male samples, one has no
real effect. For each, `benchmark/scenarios.json` states the correct decision. `python
benchmark/score.py` reports whether the pipeline reaches it.

The first run found a design flaw that no test had, because the tests encoded the same wrong
assumption as the code. An analysis with only two replicates per group passed, since a quality
check with status "warning" never blocked anything. All tests were green at the time. The fix was
a design change, in commit `73319a9`.

## In what order was it built?

`git log --reverse --oneline` shows the path. Each line is one commit: a small, complete step
with a message saying what changed.

```text
005de28 Add the payload contract skill and validator
85d5785 Add method skills for clusters, PCA, GLM, lasso, PLS-DA and synthesis
00df7c7 Add PBMC 10k worked example converted to the contract
49f8609 Add limma skill, R extractor, starter questions and simulated example
b9cc0e2 Add judge skill and Bedrock ask step with vote shares and gates
ce3f113 Add README with adaptation and HPC walkthroughs
9440a9d Add Python package with scanpy rank_genes_groups wrapper
8fd8fc8 Restructure Python into core, adapters and one shared contract file
9c19147 Add R package: core mirroring Python, limma and Seurat adapters
9181134 Add tests, fixtures, golden payloads and a CI workflow
73319a9 Add benchmark of planted problems; let each check name its blocking statuses
270913c Add technical panels and gene-set enrichment, computed in code
e357d86 Update README and skills for the packages, wrappers, panels and benchmark
ae74cc9 Add a beginner's walkthrough of how the package is built and what a release still needs
2fd04ed add licence
c7b45ed Add MIT licence and author and maintainer metadata
```

The shape to notice: definition first, then one working example on real data, then a first
wrapper, then restructuring, then tests, then a benchmark, then new features. New features came
*after* the tests, so each one could be checked against everything that already worked.

## What went wrong along the way?

Real mistakes from building this, and what caught each one:

| Mistake | What caught it |
| :--- | :--- |
| The scanpy wrapper crashed when comparing two groups, because scanpy leaves out one column in that mode | Running it on real data, not only on the case I had in mind |
| After the restructuring, the expected-features panel silently came out empty in two-group mode: an argument was not passed along | Looking at the actual output after the change. A test now asserts the panel is filled |
| A commit recorded a file deletion but none of the new files, because the staging command had failed on one bad path | Reading what the commit contained before moving on |
| Two-replicate analyses passed all gates | The benchmark |
| A new test failed although the code was right: three matching genes out of four is not statistical evidence, and the test expected a flag | Working out the number by hand (p = 0.25), then fixing the test |
| One benchmark scenario produced zero hits, so it would have been "caught" for the wrong reason | Reading the generated data before trusting the score |

Two habits come out of this list. Look at real output, not only at whether the program finished.
And when a test fails, ask which of the two is wrong, the code or the test.

## Things to try

Each takes a few minutes and shows one mechanism.

1. **Break a rule and watch a test name it.** In `python/llomics/core.py`, find `_floored` and make
   it return `0.0`. Run the tests. Read which ones fail and why. Undo it.
2. **Change a threshold in one place.** In `contract/llomics-contract.json`, change
   `min_replicates` from 3 to 5. Run `python tools/sync.py`, then the tests, and see which fail.
   The benchmark payloads are built in advance, so its score only moves after
   `Rscript benchmark/make_scenarios.R` rebuilds them. Undo the change afterwards.
3. **Add a technical panel.** Add a gene list to `contract/llomics-panels.json`, sync, and build a
   payload. It appears in both languages with no code change.
4. **Map a new tool.** Take any results table you have, rename four columns, and call
   `llomics.dge_payload()`. That is an adapter in its smallest form.

## Where do the R help pages come from?

From comments in the code. Every exported R function has a block of lines starting with `#'`
above it, describing what it does, each argument, what it returns, and an example. Running
`roxygen2::roxygenise("r/llomics")` turns those comments into the help pages in `r/llomics/man/`
and writes the `NAMESPACE` file, which lists the functions the package makes available. Neither
`man/` nor `NAMESPACE` is edited by hand.

Arguments that several functions share, such as `question` or `trust`, are described once, in
`r/llomics/R/params.R`, and each function pulls them in with `@inheritParams`. That keeps the
descriptions identical across functions. CRAN's check fails if any argument of any exported
function is left undescribed, and it caught two that slipped through on the first run.

## What is still missing for a release?

The packages work from a checkout of this repository. Publishing them so that
`pip install llomics` or `install.packages("llomics")` works needs more. Items marked *checked*
were verified on this repository; the rest are known requirements not yet attempted.

### Both packages

| Missing | Why it blocks a release |
| :--- | :--- |
| Sources for the technical gene panels | The lists in `llomics-panels.json` were compiled from general knowledge. A release needs a citation for each |
| A green CI run | The first run on GitHub passed the Python job; the R job crashed in the runner while loading a Seurat dependency. The workflow now runs the R core and the Seurat adapter as separate jobs |
| A changelog and a versioning rule | Users need to know what changed. The contract version and the two package versions are separate numbers with no stated relation |
| The name `llomics` confirmed free | Not checked on PyPI or CRAN |
| A release procedure | Tagging, building and uploading are manual and unwritten |

### Python, for PyPI

| Missing | Why it blocks a release |
| :--- | :--- |
| A successful build (*checked*: never built; this machine has no `setuptools`) | No wheel or source archive has been produced, so whether the two JSON data files are packed inside is unverified |
| An install test | The package has only been imported from the source folder. A package that works there can still fail once installed |
| Metadata in `pyproject.toml` (*checked*: licence and author are set; no `readme`, `classifiers`, `urls`) | PyPI shows an empty page and tools cannot tell which Python versions are supported |
| A README inside `python/` | The archive is built from that folder. The licence is copied in by `tools/sync.py`; the README still sits one level above |
| One source for the version (*checked*: `0.1.0` is written in both `pyproject.toml` and `__init__.py`) | Two copies will disagree after the first forgotten edit |
| The ask step inside the package | `judge_bedrock.py` and the benchmark scorer are loose scripts. An installed package would not include them, and there are no command-line entry points |
| Tests on more than one setup | Run so far on Python 3.13 and 3.11, on Windows only |
| Tested version ranges for the optional dependencies | `scanpy>=1.10` is declared; only 1.11.5 was run |
| A linter and type checker | The code has type hints, but nothing checks them |

### R, for CRAN

`R CMD check --as-cran` has been run at each step. Before the licence and author were added it
stopped at the first step. With them in place it finished with 2 warnings and 1 note: no help
pages, and an undeclared `SeuratObject`. Both warnings are now fixed, and the check ends with a
single note, `New submission`, which every first CRAN submission gets. The examples, including
the slower `\donttest` ones, and the package's own tests run inside the check.

| Missing | Why it blocks a release |
| :--- | :--- |
| `SystemRequirements` for Python | `llomics_validate()` runs a Python script, which must be declared, and the function must fail gracefully without it |
| The PDF manual | The check was run with `--no-manual` |
| Checks on other platforms | `R CMD check` has run on Windows with R 4.5.3 only. CI runs the tests, but not the check, on Linux with R 4.6. CRAN also tests macOS and the development version of R |
| A vignette and a `NEWS` file | Not strictly required, but expected of a new submission |

One point specific to this design: the R package suggests `limma`, which is distributed by
Bioconductor and not by CRAN. CRAN allows that for suggested packages, but Bioconductor may be the
more natural home for a package whose main inputs are limma fits and Seurat objects. That choice
has not been made.

## Words used here

| Word | Meaning |
| :--- | :--- |
| Payload | The small JSON file llomics produces |
| Unit | One thing a decision is made about: a cluster, or a contrast between two groups |
| Contract | The written rules a payload must follow, kept in `contract/llomics-contract.json` |
| Core | Code that builds payloads from plain tables, with no knowledge of any tool |
| Adapter | Short code that connects one tool (scanpy, Seurat, limma) to the core |
| Fixture | Input data created for tests |
| Golden file | A saved, reviewed output that tests compare against |
| Parity | R and Python producing the same result from the same input |
| Gate | A rule in code that sends a result to a person when a condition fails |
