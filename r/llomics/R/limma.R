# limma adapter: map topTable columns onto the standard table and call the core.

#' Build an llomics payload from a limma fit
#'
#' Call it after `limma::eBayes()` or `limma::treat()`. Each coefficient in `coefs` becomes one
#' contrast unit. The quality checks are computed from the fit and the sample table where possible:
#' replicates per group, the shape of the p-value histogram, group-by-batch confounding, and
#' `eBayes(trend = TRUE)`. What cannot be read from the fit, such as whether batch is in the design
#' or whether voom was used, is declared through the arguments.
#'
#' @inheritParams llomics-params
#' @inheritParams llomics_dge_payload
#' @param fit A limma `MArrayLM` object with moderated statistics.
#' @param coefs Names of the coefficients (contrasts) to report.
#' @param targets Sample table with one row per sample in the design.
#' @param group Column of `targets` holding the group of each sample.
#' @param batch Column of `targets` holding the batch of each sample. `NULL` records the batch
#'   check as `"not_run"`.
#' @param symbol_col Column of `fit$genes` holding readable feature names.
#' @param gene_sets Gene sets for enrichment of the hits: a named list of character vectors, or the
#'   path to a GMT file.
#' @param out_dir Folder for the payload (`<analysis_id>.llomics.json`) and one full table per
#'   contrast. `NULL` (the default) writes nothing.
#' @return The payload, invisibly. See [llomics_dge_payload()].
#' @export
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   set.seed(1)
#'   targets <- data.frame(group = rep(c("control", "treated"), each = 4),
#'                         batch = rep(c("b1", "b2"), times = 4))
#'   expr <- matrix(rnorm(300 * 8, mean = 7, sd = 0.4), nrow = 300,
#'                  dimnames = list(sprintf("G%03d", 1:300), NULL))
#'   expr[1:20, 5:8] <- expr[1:20, 5:8] + 2
#'   design <- stats::model.matrix(~ 0 + group + batch, data = targets)
#'   colnames(design) <- c("control", "treated", "batch_b2")
#'   fit <- limma::lmFit(expr, design)
#'   contrasts <- limma::makeContrasts(treated_vs_control = treated - control, levels = design)
#'   fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts), trend = TRUE)
#'   payload <- llomics_limma(
#'     fit, coefs = "treated_vs_control",
#'     definitions = c(treated_vs_control = "treated - control, adjusted for batch"),
#'     targets = targets, group = "group", batch = "batch", batch_in_design = TRUE,
#'     analysis_id = "example", question = "Which genes respond to treatment?",
#'     steps = c("lmFit ~ 0 + group + batch", "eBayes(trend = TRUE)"), n_features = 5)
#'   vapply(payload$trust, function(x) x$status, "")
#' }
llomics_limma <- function(fit, coefs, definitions, targets, group, analysis_id, question, steps,
                          batch = NULL, batch_in_design = NA, single_batch = FALSE, n_features = NULL,
                          panel = character(0), symbol_col = NULL, mean_variance = NULL, fdr = NULL,
                          min_effect = 0, effect_cutoff = NULL, exclude = NULL, filters = list(),
                          gene_sets = NULL, gene_sets_name = NULL, technical_panels = TRUE,
                          trust = NULL, out_dir = NULL) {
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
