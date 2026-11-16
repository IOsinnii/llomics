# llomics_limma.R
#
# Build an llomics payload (contract v0.1) from a limma fit. One unit per contrast.
# Requires the packages limma and jsonlite.
#
#   source("skills/llomics-limma/scripts/llomics_limma.R")
#   llomics_limma(
#     fit, coefs = "stim_vs_ctrl",
#     definitions = c(stim_vs_ctrl = "stimulated - control, adjusted for batch; positive logFC = higher in stimulated"),
#     targets = targets, group = "group", batch = "batch", batch_in_design = TRUE,
#     analysis_id = "stim-vs-ctrl", question = "Which genes respond to stimulation?",
#     steps = c("filterByExpr", "TMM", "voom", "lmFit ~0 + group + batch", "eBayes"),
#     mean_variance = "voom", panel = c("ISG15", "ACTB"), symbol_col = "symbol",
#     out_dir = "outputs"
#   )
#
# Writes <out_dir>/<analysis_id>.llomics.json and one full table per contrast
# (<analysis_id>.<coef>.tsv), and returns the payload invisibly.

llomics_limma <- function(fit, coefs, definitions, targets, group,
                          analysis_id, question, steps,
                          batch = NULL, batch_in_design = NA, single_batch = FALSE,
                          filters = list(), panel = character(0), symbol_col = NULL,
                          mean_variance = NULL, min_replicates = 3,
                          fdr = 0.05, lfc = 1, top_n = 15, p_floor = 1e-300,
                          out_dir = ".", write = TRUE) {
  stopifnot(requireNamespace("limma", quietly = TRUE), requireNamespace("jsonlite", quietly = TRUE))
  if (!inherits(fit, "MArrayLM")) stop("`fit` must be a limma MArrayLM object")
  if (is.null(fit$t) || is.null(fit$p.value)) stop("`fit` has no moderated statistics: run eBayes() or treat() first")
  missing_coefs <- setdiff(coefs, colnames(fit$coefficients))
  if (length(missing_coefs)) stop("coefficients not in fit: ", paste(missing_coefs, collapse = ", "))
  missing_defs <- setdiff(coefs, names(definitions))
  if (length(missing_defs)) {
    stop("`definitions` needs a named entry in words, including the direction, for: ",
         paste(missing_defs, collapse = ", "))
  }
  if (nrow(targets) != nrow(fit$design)) stop("`targets` must have one row per sample in the design")
  if (!group %in% colnames(targets)) stop("`group` is not a column of `targets`")
  if (!is.null(batch) && !batch %in% colnames(targets)) stop("`batch` is not a column of `targets`")
  if (!is.null(symbol_col) && !symbol_col %in% colnames(fit$genes)) stop("`symbol_col` is not a column of fit$genes")

  drop_null <- function(x) x[!vapply(x, is.null, logical(1))]
  check <- function(name, status, value = NULL, note = NULL) {
    drop_null(list(check = name, status = status, value = value, note = note))
  }
  worst <- function(statuses) {
    order <- c("pass", "warn", "fail")
    order[max(match(statuses, order))]
  }

  # Shape of the raw p-value histogram (10 equal bins). Heuristic thresholds, stated in the note.
  p_shape <- function(h) {
    level <- stats::median(h[5:9])
    if (max(h[2:10]) > 1.25 * h[1]) {
      return(list(status = "warn", note = sprintf(
        "bin %d holds %d p-values, more than 1.25x the first bin (%d): p-values are not concentrated near 0, which points to a misspecified model or variance",
        which.max(h[2:10]) + 1L, max(h[2:10]), h[1])))
    }
    if (h[10] > 1.5 * level) {
      return(list(status = "warn", note = sprintf(
        "excess of p-values near 1: last bin %d versus a mid-range level of %.0f", h[10], level)))
    }
    list(status = "pass", note = "non-increasing from the first bin, flat on the right")
  }

  n_floored <- 0L
  feature_rows <- function(tt, idx) {
    d <- tt[idx, , drop = FALSE]
    d <- utils::head(d[order(-abs(d$t)), , drop = FALSE], top_n)
    n_floored <<- n_floored + sum(d$adj.P.Val < p_floor)
    out <- data.frame(feature = rownames(d), stringsAsFactors = FALSE)
    if (!is.null(symbol_col)) out$symbol <- as.character(d[[symbol_col]])
    out$effect <- signif(d$logFC, 4)
    out$ave_expr <- signif(d$AveExpr, 4)
    out$stat <- signif(d$t, 4)
    out$p_adj <- signif(pmax(d$adj.P.Val, p_floor), 4)
    out
  }

  panel_state <- function(tt) {
    if (!length(panel)) return(structure(list(), names = character(0)))
    out <- lapply(panel, function(id) {
      i <- match(id, rownames(tt))
      if (is.na(i) && !is.null(symbol_col)) i <- match(id, tt[[symbol_col]])
      if (is.na(i)) return(NA)  # not measured or filtered out: null, never 0
      list(effect = signif(tt$logFC[i], 4), ave_expr = signif(tt$AveExpr[i], 4),
           stat = signif(tt$t[i], 4), p_adj = signif(max(tt$adj.P.Val[i], p_floor), 4))
    })
    stats::setNames(out, panel)
  }

  group_n <- table(as.character(targets[[group]]))
  if (write) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  units <- list()
  source_tables <- list()
  shape_statuses <- character(0)
  shape_notes <- character(0)

  for (cf in coefs) {
    tt <- limma::topTable(fit, coef = cf, number = Inf, sort.by = "none", adjust.method = "BH")
    tt <- tt[!is.na(tt$P.Value), , drop = FALSE]

    sig <- tt$adj.P.Val < fdr
    up <- sig & tt$logFC > 0
    down <- sig & tt$logFC < 0
    h <- as.integer(table(cut(tt$P.Value, breaks = seq(0, 1, 0.1), include.lowest = TRUE)))
    shape <- p_shape(h)
    shape_statuses <- c(shape_statuses, shape$status)
    shape_notes <- c(shape_notes, sprintf("%s: %s", cf, shape$note))
    pi0 <- tryCatch(limma::propTrueNull(tt$P.Value), error = function(e) NA_real_)
    q <- stats::quantile(tt$logFC, c(0.01, 0.25, 0.5, 0.75, 0.99), names = FALSE)

    rows_up <- feature_rows(tt, up)
    rows_down <- feature_rows(tt, down)

    weights <- NULL
    if (!is.null(fit$contrasts) && cf %in% colnames(fit$contrasts)) {
      w <- fit$contrasts[, cf]
      weights <- as.list(w[w != 0])
    }

    table_path <- file.path(out_dir, sprintf("%s.%s.tsv", analysis_id, make.names(cf)))
    if (write) {
      utils::write.table(cbind(feature = rownames(tt), tt), table_path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    source_tables[[cf]] <- table_path

    units[[length(units) + 1L]] <- drop_null(list(
      unit_id = cf,
      unit_type = "contrast",
      definition = unname(definitions[[cf]]),
      contrast_weights = weights,
      size = list(n_tested = nrow(tt), n_per_group = as.list(group_n)),
      facts = list(
        fdr = fdr,
        lfc_cutoff = lfc,
        n_hits = sum(sig),
        n_up = sum(up),
        n_down = sum(down),
        n_up_lfc = sum(up & tt$logFC >= lfc),
        n_down_lfc = sum(down & tt$logFC <= -lfc),
        n_near_miss = sum(tt$adj.P.Val >= fdr & tt$adj.P.Val < 2 * fdr),
        pi0 = signif(pi0, 3),
        p_histogram = h,
        logfc_quantiles = list(q01 = signif(q[1], 3), q25 = signif(q[2], 3), q50 = signif(q[3], 3),
                               q75 = signif(q[4], 3), q99 = signif(q[5], 3))
      ),
      top_features = list(up = rows_up, down = rows_down),
      panel = panel_state(tt),
      tail = list(
        n_shown = nrow(rows_up) + nrow(rows_down),
        n_eligible = sum(sig),
        rule = sprintf("adj.P.Val < %s, ranked by |t|, up to %d per direction", format(fdr), top_n)
      ),
      trust = list(check("pvalue_distribution", shape$status, note = shape$note))
    ))
  }

  # Run-level trust checks ------------------------------------------------------------------
  min_n <- min(group_n)
  rep_status <- if (min_n >= min_replicates) "pass" else if (min_n >= 2) "warn" else "fail"
  rep_note <- sprintf("samples per group: %s; pass needs at least %d",
                      paste(sprintf("%s=%d", names(group_n), as.integer(group_n)), collapse = ", "), min_replicates)

  if (single_batch) {
    batch_check <- check("batch_confounding", "not_applicable", note = "declared by the analyst as a single batch")
  } else if (is.null(batch)) {
    batch_check <- check("batch_confounding", "not_run", note = "no batch variable was supplied to the extractor")
  } else {
    ct <- table(as.character(targets[[group]]), as.character(targets[[batch]]))
    crosstab <- stats::setNames(lapply(rownames(ct), function(g) {
      as.list(stats::setNames(as.integer(ct[g, ]), colnames(ct)))
    }), rownames(ct))
    if (all(colSums(ct > 0) == 1)) {
      batch_check <- check("batch_confounding", "fail", crosstab,
                           "every batch contains a single group: the group effect cannot be separated from batch")
    } else if (!isTRUE(batch_in_design)) {
      batch_check <- check("batch_confounding", "warn", crosstab,
                           if (is.na(batch_in_design)) "not declared whether batch is in the design (set batch_in_design)"
                           else "batch is not in the design")
    } else if (any(ct == 0)) {
      batch_check <- check("batch_confounding", "warn", crosstab,
                           "batch is in the design, but some group x batch cells are empty")
    } else {
      batch_check <- check("batch_confounding", "pass", crosstab, "batch is in the design and every group occurs in every batch")
    }
  }

  mt_check <- if (length(coefs) == 1) {
    check("multiple_testing", "pass", note = "Benjamini-Hochberg within the one contrast reported")
  } else {
    check("multiple_testing", "warn", value = length(coefs), note = sprintf(
      "Benjamini-Hochberg within each of %d contrasts; not adjusted across contrasts", length(coefs)))
  }

  trend_detected <- length(fit$s2.prior) > 1
  mv_check <- if (!is.null(mean_variance)) {
    switch(mean_variance,
      none = check("mean_variance", "warn", note = "declared: no mean-variance modelling"),
      not_applicable = check("mean_variance", "not_applicable", note = "declared not applicable to this data type"),
      check("mean_variance", "pass", note = paste("declared:", mean_variance)))
  } else if (trend_detected) {
    check("mean_variance", "pass", note = "eBayes(trend = TRUE) detected from the fit")
  } else {
    check("mean_variance", "not_run", note = "not declared; pass mean_variance = 'voom', 'trend', 'none' or 'not_applicable'")
  }

  df_prior <- fit$df.prior
  payload <- list(
    llomics = "0.1",
    method = "limma",
    run = list(
      analysis_id = analysis_id,
      question = question,
      tool = paste("limma", as.character(utils::packageVersion("limma"))),
      steps = I(as.character(steps)),
      design_columns = I(colnames(fit$design)),
      filters = c(filters, list(
        list(rule = sprintf("adj.P.Val < %s (Benjamini-Hochberg)", format(fdr)), why = "FDR control"),
        list(rule = sprintf("top %d per direction by |t|; never padded", top_n), why = "payload size")
      )),
      p_floor = p_floor,
      n_p_values_floored = n_floored,
      df_prior_median = if (all(is.finite(df_prior))) signif(stats::median(df_prior), 4) else NA,
      df_prior_infinite = any(is.infinite(df_prior)),
      source_tables = source_tables
    ),
    dataset = list(
      n_observations = nrow(fit$design),
      observation_type = "sample",
      n_features = nrow(fit$coefficients),
      reference = list(n_contrasts = length(coefs))
    ),
    trust = list(
      check("replicates", rep_status, value = as.integer(min_n), note = rep_note),
      check("pvalue_distribution", worst(shape_statuses), note = paste(shape_notes, collapse = " | ")),
      batch_check,
      mt_check,
      mv_check
    ),
    units = units
  )

  if (write) {
    path <- file.path(out_dir, paste0(analysis_id, ".llomics.json"))
    jsonlite::write_json(payload, path, auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE)
    message("Wrote ", path, " (", length(units), " unit", if (length(units) != 1) "s", ")")
  }
  invisible(payload)
}
