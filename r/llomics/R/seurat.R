# Seurat adapters: wrap FindAllMarkers and FindMarkers, map their columns onto the standard table,
# and call the core.

.seurat_tool <- function(fun, dots) {
  args <- if (length(dots)) {
    paste(sprintf("%s = %s", names(dots), vapply(dots, function(v) paste(deparse(v), collapse = ""), character(1))),
          collapse = ", ")
  } else {
    "defaults"
  }
  sprintf("Seurat %s %s(%s)", as.character(utils::packageVersion("Seurat")), fun, args)
}

.seurat_steps <- function(object, steps) {
  if (!is.null(steps)) return(list(steps = steps, source = "supplied by the analyst"))
  list(steps = names(object@commands),
       source = "names of the commands logged in the Seurat object; cell and gene QC filters are not recorded there")
}

.as_doublet <- function(x) if (is.logical(x)) x else tolower(as.character(x)) %in% c("true", "doublet")

.meta_col <- function(meta, col, what) {
  if (is.null(col)) return(NULL)
  if (!col %in% colnames(meta)) stop(what, " '", col, "' is not a metadata column of the object")
  meta[[col]]
}

.seurat_prefilter <- function(dots) {
  # Seurat 5 defaults: logfc.threshold = 0.1, min.pct = 0.01
  lfc <- dots$logfc.threshold %||% 0.1
  pct <- dots$min.pct %||% 0.01
  list(active = lfc > 0 || pct > 0,
       note = sprintf("Seurat tested only features passing logfc.threshold = %s and min.pct = %s", .fmt(lfc), .fmt(pct)))
}

#' Run Seurat::FindAllMarkers and build an llomics payload from the result
#'
#' Each identity class becomes a cluster unit with its top markers, the fraction of cells
#' expressing each one inside and outside the cluster, its quality-control medians compared with
#' the whole dataset, its nearest cluster in PCA space, and the mean expression of the panel
#' features. Seurat reports no test statistic, so markers are ranked by adjusted p-value and then
#' fold change.
#'
#' @inheritParams llomics-params
#' @param object A Seurat object with identities set.
#' @param n_features Markers loaded per cluster. Never padded. `NULL` uses the contract default (15).
#' @param panel Features whose mean expression is reported for every cluster. A feature that is not
#'   in the object is written as `null`.
#' @param analysis_id Short identifier for the analysis, stored in the run card.
#' @param out Path for the payload JSON. The full marker table is written next to it. `NULL`
#'   writes nothing.
#' @param min_logfc Log fold change a marker must exceed to be eligible for the top list.
#' @param exclude Regular expression. Markers whose names match are left out of the top list. The
#'   default excludes mitochondrial, ribosomal-protein and MALAT1 genes.
#' @param sample_col,batch_col,doublet_col Metadata columns holding the biological sample, the
#'   batch and the predicted-doublet flag. They feed the quality checks; a missing one records the
#'   related check as `"not_run"`.
#' @param qc_cols Metadata columns summarised per cluster and compared with the dataset median.
#'   Columns that are absent are skipped.
#' @param steps Character vector describing the analysis. `NULL` uses the names of the commands
#'   logged in the Seurat object.
#' @param markers An existing `FindAllMarkers()` result to use instead of running the test.
#' @param write_table Also write the full marker table next to `out`.
#' @param gene_sets Gene sets for enrichment of the markers: a named list of character vectors, or
#'   the path to a GMT file.
#' @param ... Passed to `Seurat::FindAllMarkers()`.
#' @return The payload, invisibly. See [llomics_clusters_payload()].
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("Seurat", quietly = TRUE) &&
#'     requireNamespace("SeuratObject", quietly = TRUE)) {
#'   data("pbmc_small", package = "SeuratObject")
#'   payload <- llomics_FindAllMarkers(pbmc_small, n_features = 5,
#'                                     panel = c("CD3D", "MS4A1", "LYZ"),
#'                                     question = "Which cell type is each cluster?",
#'                                     verbose = FALSE)
#'   payload$units[[1]]$top_features
#' }
#' }
llomics_FindAllMarkers <- function(object, n_features = NULL, panel = character(0), question = NULL,
                                   analysis_id = NULL, out = NULL, fdr = NULL, min_logfc = 0,
                                   exclude = llomics_contract()$defaults$exclude, sample_col = NULL,
                                   batch_col = NULL, doublet_col = NULL,
                                   qc_cols = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                   steps = NULL, trust = NULL, markers = NULL, write_table = TRUE,
                                   gene_sets = NULL, gene_sets_name = NULL, technical_panels = TRUE, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("the Seurat package is needed for llomics_FindAllMarkers()")
  dots <- list(...)
  if (is.character(gene_sets) && length(gene_sets) == 1) {
    gene_sets_name <- gene_sets_name %||% basename(gene_sets)
    gene_sets <- llomics_read_gmt(gene_sets)
  }
  if (is.null(markers)) markers <- Seurat::FindAllMarkers(object, ...)
  if (!nrow(markers)) stop("FindAllMarkers returned no rows")

  idents <- Seurat::Idents(object)
  clusters <- levels(idents)
  labels <- as.character(idents)
  meta <- object[[]]

  tables <- stats::setNames(lapply(clusters, function(cl) {
    m <- markers[as.character(markers$cluster) == cl, , drop = FALSE]
    data.frame(feature = as.character(m$gene), effect = m$avg_log2FC, p = m$p_val, p_adj = m$p_val_adj,
               pct_in = m$pct.1, pct_out = m$pct.2, stringsAsFactors = FALSE)
  }), clusters)

  present <- intersect(panel, rownames(object))
  panel_values <- stats::setNames(lapply(clusters, function(cl) {
    stats::setNames(as.list(rep(NA, length(panel))), panel)
  }), clusters)
  if (length(present)) {
    expr <- Seurat::FetchData(object, vars = present)
    for (cl in clusters) {
      for (g in present) panel_values[[cl]][[g]] <- .sig(mean(expr[labels == cl, g]))
    }
  }

  extra_facts <- NULL
  if ("pca" %in% names(object@reductions) && length(clusters) > 2) {
    emb <- Seurat::Embeddings(object, "pca")
    cent <- rowsum(emb, labels) / as.vector(table(labels)[rownames(rowsum(emb, labels))])
    corr <- stats::cor(t(cent))
    diag(corr) <- -Inf
    extra_facts <- stats::setNames(lapply(rownames(corr), function(cl) {
      list(nearest_cluster = colnames(corr)[which.max(corr[cl, ])])
    }), rownames(corr))
  }

  source_tables <- list(object = "FindAllMarkers result passed to llomics")
  if (!is.null(out) && write_table) {
    table_path <- paste0(sub("\\.json$", "", out), ".table.tsv")
    dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(markers, table_path, sep = "\t", quote = FALSE, row.names = FALSE)
    source_tables$table <- table_path
  }

  st <- .seurat_steps(object, steps)
  pre <- .seurat_prefilter(dots)
  qc_present <- intersect(qc_cols, colnames(meta))
  doublet <- .meta_col(meta, doublet_col, "doublet_col")
  payload <- llomics_clusters_payload(
    tables,
    analysis_id = analysis_id %||% "seurat-FindAllMarkers", question = question,
    tool = .seurat_tool("FindAllMarkers", dots), steps = st$steps,
    test = paste(dots$test.use %||% "wilcox", "test"),
    labels = labels,
    qc = if (length(qc_present)) as.list(meta[qc_present]),
    doublet = if (!is.null(doublet)) .as_doublet(doublet),
    batch = .meta_col(meta, batch_col, "batch_col"),
    sample = .meta_col(meta, sample_col, "sample_col"),
    panel_values = panel_values, extra_facts = extra_facts, n_dataset_features = nrow(object),
    n_features = n_features, fdr = fdr, min_effect = min_logfc, exclude = exclude,
    filters = if (pre$active) list(list(rule = pre$note, why = "Seurat pre-filter; n_eligible counts features among those tested")) else list(),
    source_tables = source_tables,
    run_extra = list(steps_source = st$source,
                     p_adj_note = "Seurat p_val_adj is a Bonferroni correction over all features in the object"),
    gene_sets = gene_sets, gene_sets_name = gene_sets_name, technical_panels = technical_panels,
    trust = trust
  )
  if (!is.null(out)) llomics_write(payload, out)
  invisible(payload)
}

#' Run Seurat::FindMarkers and build an llomics payload from the result
#'
#' Builds one contrast unit: `ident.1` against `ident.2`, or against all other cells. By default
#' Seurat tests only features passing `logfc.threshold` and `min.pct`; the p-value histogram is then
#' not a sample of all features and its check is recorded as `"not_run"`. Pass
#' `logfc.threshold = 0, min.pct = 0` to test every feature.
#'
#' @inheritParams llomics_FindAllMarkers
#' @param ident.1,ident.2,group.by As in `Seurat::FindMarkers()`. With `ident.2 = NULL` the
#'   comparison is against all other cells.
#' @param n_features Features loaded per direction. Never padded. `NULL` uses the contract default
#'   (15).
#' @param panel Features reported with their statistics whether or not they rank.
#' @param markers An existing `FindMarkers()` result to use instead of running the test.
#' @param sample_col,batch_col Metadata columns holding the biological sample and the batch of each
#'   cell. They feed the replicates and batch checks; a missing one records the check as `"not_run"`.
#' @param ... Passed to `Seurat::FindMarkers()`.
#' @return The payload, invisibly. See [llomics_dge_payload()].
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("Seurat", quietly = TRUE) &&
#'     requireNamespace("SeuratObject", quietly = TRUE)) {
#'   data("pbmc_small", package = "SeuratObject")
#'   payload <- llomics_FindMarkers(pbmc_small, ident.1 = "0", ident.2 = "1",
#'                                  logfc.threshold = 0, min.pct = 0, n_features = 5,
#'                                  question = "How do clusters 0 and 1 differ?",
#'                                  verbose = FALSE)
#'   payload$units[[1]]$top_features
#' }
#' }
llomics_FindMarkers <- function(object, ident.1, ident.2 = NULL, group.by = NULL, n_features = NULL,
                                panel = character(0), question = NULL, analysis_id = NULL, out = NULL,
                                fdr = NULL, min_logfc = 0, exclude = llomics_contract()$defaults$exclude,
                                sample_col = NULL, batch_col = NULL, steps = NULL, trust = NULL,
                                markers = NULL, write_table = TRUE, gene_sets = NULL, gene_sets_name = NULL,
                                technical_panels = TRUE, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("the Seurat package is needed for llomics_FindMarkers()")
  dots <- list(...)
  if (is.character(gene_sets) && length(gene_sets) == 1) {
    gene_sets_name <- gene_sets_name %||% basename(gene_sets)
    gene_sets <- llomics_read_gmt(gene_sets)
  }
  if (is.null(markers)) {
    markers <- Seurat::FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2, group.by = group.by, ...)
  }
  meta <- object[[]]
  labels <- if (!is.null(group.by)) as.character(meta[[group.by]]) else as.character(Seurat::Idents(object))
  name1 <- paste(ident.1, collapse = "+")
  name2 <- if (is.null(ident.2)) "rest" else paste(ident.2, collapse = "+")
  in1 <- labels %in% ident.1
  in2 <- if (is.null(ident.2)) !in1 else labels %in% ident.2
  uid <- paste0(name1, "_vs_", name2)

  tab <- data.frame(feature = rownames(markers), effect = markers$avg_log2FC, p = markers$p_val,
                    p_adj = markers$p_val_adj, pct_in = markers$pct.1, pct_out = markers$pct.2,
                    stringsAsFactors = FALSE)

  size <- list(n_cells = stats::setNames(list(sum(in1), sum(in2)), c(name1, name2)))
  sample <- .meta_col(meta, sample_col, "sample_col")
  batch <- .meta_col(meta, batch_col, "batch_col")
  group_sizes <- NULL
  if (!is.null(sample)) {
    group_sizes <- stats::setNames(list(length(unique(sample[in1])), length(unique(sample[in2]))), c(name1, name2))
    size$n_samples <- group_sizes
  }
  crosstab <- NULL
  if (!is.null(batch)) {
    keep <- in1 | in2
    crosstab <- .crosstab(ifelse(in1[keep], name1, name2), batch[keep])
  }
  batch_modelled <- !is.null(batch_col) && batch_col %in% dots$latent.vars

  test <- dots$test.use %||% "wilcox"
  mean_variance <- if (test %in% c("wilcox", "wilcox_limma", "roc")) "not_applicable" else
    if (test %in% c("t", "bimod")) "none" else test

  source_tables <- list(object = "FindMarkers result passed to llomics")
  if (!is.null(out) && write_table) {
    table_path <- paste0(sub("\\.json$", "", out), ".table.tsv")
    dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(cbind(feature = rownames(markers), markers), table_path, sep = "\t", quote = FALSE, row.names = FALSE)
    source_tables$table <- table_path
  }

  st <- .seurat_steps(object, steps)
  pre <- .seurat_prefilter(dots)
  where <- if (!is.null(group.by)) group.by else "the active identities"
  payload <- llomics_dge_payload(
    stats::setNames(list(tab), uid),
    stats::setNames(sprintf("%s - %s in %s; positive effect = higher in %s", name1, name2, where, name1), uid),
    analysis_id = analysis_id %||% paste0("seurat-", uid), question = question,
    tool = .seurat_tool("FindMarkers", dots), steps = st$steps,
    n_observations = ncol(object), observation_type = "cell",
    group_sizes = group_sizes, units_are_replicates = FALSE,
    batch_crosstab = crosstab, batch_in_design = if (batch_modelled) TRUE else NA,
    batch_adjustable = batch_modelled,
    mean_variance = mean_variance, correction = "Bonferroni over all features in the object",
    n_features = n_features, fdr = fdr, min_effect = min_logfc, exclude = exclude, panel = panel,
    sizes = stats::setNames(list(size), uid),
    filters = if (pre$active) list(list(rule = pre$note, why = "Seurat pre-filter; counts refer to the features tested")) else list(),
    source_tables = source_tables,
    run_extra = list(steps_source = st$source),
    dataset_reference = list(reference_group = name2),
    p_shape_skip = if (pre$active) paste0(pre$note, ", so the p-value histogram is not a sample of all features; rerun with logfc.threshold = 0 and min.pct = 0 to judge it"),
    gene_sets = gene_sets, gene_sets_name = gene_sets_name, technical_panels = technical_panels,
    trust = trust
  )
  if (!is.null(out)) llomics_write(payload, out)
  invisible(payload)
}
