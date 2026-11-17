# limma adapter: map topTable columns onto the standard table and call the core.

#' Build an llomics payload from a limma fit
#'
#' Call it after `eBayes()` (or `treat()`). One unit per coefficient in `coefs`.
#'
#' @param fit A limma `MArrayLM` with moderated statistics.
#' @param coefs Coefficient (contrast) names to report.
#' @param definitions Named character vector: each contrast in words, including the direction,
#'   e.g. `c(stim_vs_ctrl = "stimulated - control, adjusted for batch; positive effect = higher in stimulated")`.
#' @param targets Sample table, one row per sample in the design.
#' @param group,batch Column names in `targets`.
#' @param batch_in_design Is batch a term in the design? Declare it; it cannot be read from the fit.
#' @param n_features How many features to load per direction.
#' @param panel Feature IDs or symbols reported whether or not they rank.
#' @param symbol_col Column of `fit$genes` holding readable names.
#' @param mean_variance One of "voom", "trend", "none", "not_applicable". `eBayes(trend = TRUE)` is detected.
#' @param gene_sets Named list of gene sets, or the path to a GMT file, for enrichment of the hits.
#' @param technical_panels TRUE scores the shipped artefact panels against the hits; FALSE skips it.
#' @param out_dir Where the payload and the full table per contrast are written. NULL writes nothing.
#' @export
llomics_limma <- function(fit, coefs, definitions, targets, group, analysis_id, question, steps,
                          batch = NULL, batch_in_design = NA, single_batch = FALSE, n_features = NULL,
                          panel = character(0), symbol_col = NULL, mean_variance = NULL, fdr = NULL,
                          min_effect = 0, effect_cutoff = NULL, exclude = NULL, filters = list(),
                          gene_sets = NULL, gene_sets_name = NULL, technical_panels = TRUE,
                          trust = NULL, out_dir = ".") {
  if (!requireNamespace("limma", quietly = TRUE)) stop("the limma package is needed for llomics_limma()")
  if (!inherits(fit, "MArrayLM")) stop("`fit` must be a limma MArrayLM object")
  if (is.null(fit$t) || is.null(fit$p.value)) stop("`fit` has no moderated statistics: run eBayes() or treat() first")
  missing_coefs <- setdiff(coefs, colnames(fit$coefficients))
  if (length(missing_coefs)) stop("coefficients not in fit: ", paste(missing_coefs, collapse = ", "))
  if (nrow(targets) != nrow(fit$design)) stop("`targets` must have one row per sample in the design")
  for (col in c(group, batch)) if (!col %in% colnames(targets)) stop("`", col, "` is not a column of `targets`")
  if (!is.null(symbol_col) && !symbol_col %in% colnames(fit$genes)) stop("`symbol_col` is not a column of fit$genes")

  if (is.character(gene_sets) && length(gene_sets) == 1) {
    gene_sets_name <- gene_sets_name %||% basename(gene_sets)
    gene_sets <- llomics_read_gmt(gene_sets)
  }
  group_n <- table(as.character(targets[[group]]))
  group_sizes <- as.list(stats::setNames(as.integer(group_n), names(group_n)))

  tables <- list()
  sizes <- list()
  weights <- list()
  source_tables <- list()
  if (!is.null(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (cf in coefs) {
    tt <- limma::topTable(fit, coef = cf, number = Inf, sort.by = "none", adjust.method = "BH")
    tab <- data.frame(feature = rownames(tt), stringsAsFactors = FALSE)
    if (!is.null(symbol_col)) tab$symbol <- as.character(tt[[symbol_col]])
    tab$effect <- tt$logFC
    tab$ave_expr <- tt$AveExpr
    tab$stat <- tt$t
    tab$p <- tt$P.Value
    tab$p_adj <- tt$adj.P.Val
    tables[[cf]] <- tab
    sizes[[cf]] <- list(n_per_group = group_sizes)
    if (!is.null(fit$contrasts) && cf %in% colnames(fit$contrasts)) {
      w <- fit$contrasts[, cf]
      weights[[cf]] <- as.list(w[w != 0])
    }
    if (!is.null(out_dir)) {
      path <- file.path(out_dir, sprintf("%s.%s.tsv", analysis_id, make.names(cf)))
      utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
      source_tables[[cf]] <- path
    }
  }

  df_prior <- fit$df.prior
  payload <- llomics_dge_payload(
    tables, definitions,
    analysis_id = analysis_id, question = question,
    tool = paste("limma", as.character(utils::packageVersion("limma"))),
    steps = steps, n_observations = nrow(fit$design), observation_type = "sample",
    group_sizes = group_sizes, units_are_replicates = TRUE,
    batch_crosstab = if (!is.null(batch)) .crosstab(targets[[group]], targets[[batch]]),
    batch_in_design = batch_in_design, batch_adjustable = TRUE, single_batch = single_batch,
    mean_variance = mean_variance,
    mean_variance_detected = if (length(fit$s2.prior) > 1) "eBayes(trend = TRUE)",
    n_features = n_features, fdr = fdr, min_effect = min_effect, effect_cutoff = effect_cutoff,
    exclude = exclude, panel = panel, sizes = sizes, contrast_weights = weights, filters = filters,
    source_tables = source_tables,
    run_extra = list(
      design_columns = I(colnames(fit$design)),
      df_prior_median = if (all(is.finite(df_prior))) .sig(stats::median(df_prior)) else NA,
      df_prior_infinite = any(is.infinite(df_prior))
    ),
    gene_sets = gene_sets, gene_sets_name = gene_sets_name, technical_panels = technical_panels,
    trust = trust
  )
  if (!is.null(out_dir)) llomics_write(payload, file.path(out_dir, paste0(analysis_id, ".llomics.json")))
  invisible(payload)
}
