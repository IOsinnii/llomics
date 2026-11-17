# Build the fixture payloads with the R core. Called by tests/test_parity.py:
#
#   Rscript tests/parity/build_from_fixtures.R <fixtures_dir> <out_dir>
#
# Reads the same tables and argument files as tests/helpers.py and writes dge.json and clusters.json.

suppressPackageStartupMessages(library(llomics))
args <- commandArgs(trailingOnly = TRUE)
fixtures <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  utils::read.delim(path, na.strings = "NA", stringsAsFactors = FALSE, check.names = FALSE, quote = "",
                    colClasses = c(feature = "character"))
}
read_args <- function(name) jsonlite::fromJSON(file.path(fixtures, name), simplifyVector = FALSE)
chr <- function(x) unlist(x, use.names = FALSE)

# ---- dge ----
a <- read_args("dge_args.json")
uids <- names(a$definitions)
tables <- stats::setNames(lapply(uids, function(u) read_table(file.path(fixtures, sprintf("dge_%s.tsv", u)))), uids)
dge <- llomics_dge_payload(
  tables, a$definitions,
  analysis_id = a$analysis_id, question = a$question, tool = a$tool, steps = chr(a$steps),
  n_observations = a$n_observations, group_sizes = a$group_sizes, batch_crosstab = a$batch_crosstab,
  batch_in_design = a$batch_in_design, mean_variance = a$mean_variance, n_features = a$n_features,
  fdr = a$fdr, min_effect = a$min_effect, exclude = a$exclude, panel = chr(a$panel), sizes = a$sizes,
  contrast_weights = a$contrast_weights, filters = a$filters, source_tables = a$source_tables, trust = a$trust)
llomics_write(dge, file.path(out_dir, "dge.json"))

# ---- clusters ----
a <- read_args("clusters_args.json")
markers <- utils::read.delim(file.path(fixtures, "markers.tsv"), na.strings = "NA", stringsAsFactors = FALSE, quote = "")
clusters <- unique(markers$cluster)
tables <- stats::setNames(lapply(clusters, function(cl) {
  markers[markers$cluster == cl, setdiff(names(markers), "cluster"), drop = FALSE]
}), clusters)
obs <- utils::read.delim(file.path(fixtures, "obs.tsv"), stringsAsFactors = FALSE, quote = "")
panel_values <- lapply(a$panel_values, function(p) lapply(p, function(v) if (is.null(v)) NA else v))
cl <- llomics_clusters_payload(
  tables,
  analysis_id = a$analysis_id, question = a$question, tool = a$tool, steps = chr(a$steps), test = a$test,
  labels = obs$label, qc = list(n_genes = obs$n_genes, pct_mt = obs$pct_mt), doublet = obs$doublet,
  batch = obs$batch, sample = obs$sample, panel_values = panel_values, extra_facts = a$extra_facts,
  n_dataset_features = a$n_dataset_features, n_features = a$n_features, fdr = a$fdr,
  min_effect = a$min_effect, exclude = a$exclude, source_tables = a$source_tables)
llomics_write(cl, file.path(out_dir, "clusters.json"))
