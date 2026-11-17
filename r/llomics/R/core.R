# Core of llomics: plain tables in, payload units out. Mirrors python/llomics/core.py.
#
# A table is a data.frame. Required columns: feature, effect, p, p_adj. Optional: symbol, stat,
# ave_expr, pct_in, pct_out. Adapters for limma, Seurat and others only map their columns onto
# these names and collect provenance; everything that decides what goes into a payload lives here.
# Notes and rules are worded exactly as in the Python core, and tests/parity checks that the two
# produce the same payload from the same table.

.contract_env <- new.env(parent = emptyenv())

#' The llomics contract shared by R, Python and the validator
#' @export
llomics_contract <- function() {
  if (is.null(.contract_env$contract)) {
    path <- system.file("llomics-contract.json", package = "llomics")
    if (!nzchar(path)) stop("llomics-contract.json not found: run python tools/sync.py, then reinstall the package")
    .contract_env$contract <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  }
  .contract_env$contract
}

.default <- function(name) llomics_contract()$defaults[[name]]
`%||%` <- function(a, b) if (is.null(a)) b else a
.fmt <- function(x) format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
.empty_object <- function() structure(list(), names = character(0))
.drop_null <- function(x) x[!vapply(x, is.null, logical(1))]

.sig <- function(x, digits = NULL) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  signif(x, digits %||% .default("digits"))
}

# ---- trust checks ----------------------------------------------------------------------------

.check <- function(name, status, note = NULL, value = NULL) {
  .drop_null(list(check = name, status = status, value = value, note = note))
}

.worst <- function(statuses) {
  ord <- c("pass", "warn", "fail")
  ord[max(match(statuses, ord))]
}

.apply_trust_overrides <- function(checks, overrides) {
  if (!length(overrides)) return(checks)
  names(checks) <- vapply(checks, function(c) c$check, character(1))
  for (name in names(overrides)) {
    decl <- overrides[[name]]
    status <- if (is.list(decl)) decl$status else decl[[1]]
    note <- if (is.list(decl)) decl$note else if (length(decl) > 1) decl[[2]] else NULL
    checks[[name]] <- .check(name, status,
                             if (!is.null(note)) paste0("declared by the analyst: ", note) else "declared by the analyst")
  }
  unname(checks)
}

.p_histogram <- function(p) {
  p <- p[is.finite(p)]
  tabulate(pmin(floor(p * 10), 9) + 1, nbins = 10)
}

.p_shape <- function(h) {
  cfg <- llomics_contract()$p_shape
  if (sum(h) == 0) return(list(status = "warn", note = "no p-values"))
  mid <- sort(h[5:9])[3]
  later <- max(h[2:10])
  if (later > cfg$later_bin_ratio * h[1]) {
    return(list(status = "warn", note = sprintf(
      "bin %d holds %d p-values, more than %sx the first bin (%d): p-values are not concentrated near 0, which points to a misspecified test or variance",
      which.max(h[2:10]) + 1L, later, .fmt(cfg$later_bin_ratio), h[1])))
  }
  if (h[10] > cfg$last_bin_ratio * mid) {
    return(list(status = "warn", note = sprintf(
      "excess of p-values near 1: last bin %d versus a mid-range level of %d", h[10], mid)))
  }
  list(status = "pass", note = "non-increasing from the first bin, flat on the right")
}

.storey_pi0 <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  min(1, mean(p > 0.5) / 0.5)
}

.replicates_check <- function(group_sizes, min_replicates = NULL, units_are_replicates = TRUE, what = "samples") {
  min_replicates <- min_replicates %||% .default("min_replicates")
  if (!length(group_sizes)) {
    note <- sprintf("number of %s per group not supplied", what)
    if (!units_are_replicates) {
      note <- paste0(note, "; the test units were not independent replicates, so p-values do not reflect variation between samples")
    }
    return(.check("replicates", "not_run", note))
  }
  smallest <- min(unlist(group_sizes))
  listing <- paste(sprintf("%s=%d", names(group_sizes), as.integer(unlist(group_sizes))), collapse = ", ")
  if (!units_are_replicates) {
    status <- if (smallest < 2) "fail" else "warn"
    note <- sprintf("%s per group: %s. The test units were not these replicates, so p-values do not reflect variation between samples; aggregate per sample (pseudobulk) for inference", what, listing)
  } else {
    status <- if (smallest >= min_replicates) "pass" else if (smallest >= 2) "warn" else "fail"
    note <- sprintf("%s per group: %s; pass needs at least %d", what, listing, as.integer(min_replicates))
  }
  .check("replicates", status, note, as.integer(smallest))
}

# crosstab: named list by group of named lists by batch, holding observation counts
.batch_confounding_check <- function(crosstab, batch_in_design = NA, adjustable = TRUE, single_batch = FALSE) {
  name <- "batch_confounding"
  if (single_batch) return(.check(name, "not_applicable", "declared by the analyst as a single batch"))
  if (!length(crosstab)) return(.check(name, "not_run", "no batch variable was supplied"))
  count <- function(row, b) if (is.null(row[[b]])) 0 else row[[b]]
  all_batches <- unique(unlist(lapply(crosstab, names)))
  batches <- sort(all_batches[vapply(all_batches, function(b) any(vapply(crosstab, count, numeric(1), b = b) > 0), logical(1))],
                  method = "radix")
  if (length(batches) <= 1) return(.check(name, "not_applicable", "all observations come from one batch"))
  groups_per_batch <- vapply(batches, function(b) sum(vapply(crosstab, count, numeric(1), b = b) > 0), numeric(1))
  if (all(groups_per_batch == 1)) {
    return(.check(name, "fail", "every batch contains a single group: the group effect cannot be separated from batch", crosstab))
  }
  if (!adjustable) return(.check(name, "warn", "the test has no covariates, so batch is not adjusted", crosstab))
  if (!isTRUE(batch_in_design)) {
    note <- if (is.na(batch_in_design)) "not declared whether batch is in the design (set batch_in_design)" else "batch is not in the design"
    return(.check(name, "warn", note, crosstab))
  }
  empty <- any(vapply(crosstab, function(row) any(vapply(batches, function(b) count(row, b) == 0, logical(1))), logical(1)))
  if (empty) return(.check(name, "warn", "batch is in the design, but some group x batch cells are empty", crosstab))
  .check(name, "pass", "batch is in the design and every group occurs in every batch", crosstab)
}

.crosstab <- function(groups, batches) {
  ct <- table(as.character(groups), as.character(batches))
  stats::setNames(lapply(rownames(ct), function(g) {
    row <- ct[g, ]
    as.list(stats::setNames(as.integer(row), colnames(ct)))[row > 0]
  }), rownames(ct))
}

.multiple_testing_check <- function(n_contrasts, method = "Benjamini-Hochberg") {
  if (n_contrasts == 1) return(.check("multiple_testing", "pass", sprintf("%s within the one contrast reported", method)))
  .check("multiple_testing", "warn",
         sprintf("%s within each of %d contrasts; not adjusted across contrasts", method, as.integer(n_contrasts)),
         as.integer(n_contrasts))
}

.mean_variance_check <- function(declared = NULL, detected = NULL) {
  name <- "mean_variance"
  if (identical(declared, "none")) return(.check(name, "warn", "declared: no mean-variance modelling"))
  if (identical(declared, "not_applicable")) return(.check(name, "not_applicable", "declared not applicable to this test or data type"))
  if (!is.null(declared)) return(.check(name, "pass", paste0("declared: ", declared)))
  if (!is.null(detected)) return(.check(name, "pass", paste0("detected: ", detected)))
  .check(name, "not_run", "not declared; set mean_variance to 'voom', 'trend', 'none' or 'not_applicable'")
}

# ---- units -----------------------------------------------------------------------------------

.rows <- function(table) {
  required <- llomics_contract()$feature_columns$required
  missing <- setdiff(required, names(table))
  if (length(missing)) stop("table is missing the required column(s): ", paste(missing, collapse = ", "))
  table <- as.data.frame(table, stringsAsFactors = FALSE)
  table$feature <- as.character(table$feature)
  table <- table[is.finite(table$p) & is.finite(table$effect), , drop = FALSE]
  table$p_adj[!is.finite(table$p_adj)] <- NA_real_
  if ("stat" %in% names(table)) table$stat[!is.finite(table$stat)] <- 0
  if ("symbol" %in% names(table)) {
    table$symbol <- as.character(table$symbol)
    table$symbol[!is.na(table$symbol) & table$symbol == ""] <- NA_character_
  }
  rownames(table) <- NULL
  table
}

.label <- function(rows) {
  if ("symbol" %in% names(rows)) ifelse(is.na(rows$symbol), rows$feature, rows$symbol) else rows$feature
}

.excluded <- function(rows, exclude) {
  if (is.null(exclude) || !nrow(rows)) return(rep(FALSE, nrow(rows)))
  grepl(exclude, .label(rows), perl = TRUE)
}

.below <- function(p_adj, cutoff) !is.na(p_adj) & p_adj < cutoff

.feature_rows <- function(rows, p_floor) {
  out <- data.frame(feature = rows$feature, stringsAsFactors = FALSE)
  if ("symbol" %in% names(rows)) out$symbol <- rows$symbol
  out$effect <- .sig(rows$effect)
  if ("ave_expr" %in% names(rows)) out$ave_expr <- .sig(rows$ave_expr)
  if ("stat" %in% names(rows)) out$stat <- .sig(rows$stat)
  out$p_adj <- .sig(pmax(rows$p_adj, p_floor))
  if ("pct_in" %in% names(rows)) out$pct_in <- .sig(rows$pct_in, 3)
  if ("pct_out" %in% names(rows)) out$pct_out <- .sig(rows$pct_out, 3)
  out
}

# Rank by test statistic where there is one, else by p-value then effect size. Stable for ties.
.rank <- function(rows, has_stat, signed = 0) {
  if (!nrow(rows)) return(rows)
  o <- if (has_stat) {
    if (signed == 0) order(-abs(rows$stat), method = "radix") else order(-signed * rows$stat, method = "radix")
  } else {
    order(rows$p, -abs(rows$effect), method = "radix")
  }
  rows[o, , drop = FALSE]
}

.count_floored <- function(table, p_floor) sum(is.finite(table$p_adj) & table$p_adj < p_floor)

.contrast_unit <- function(table, unit_id, definition, n_features = NULL, fdr = NULL, min_effect = 0,
                           effect_cutoff = NULL, exclude = NULL, panel = character(0), size = NULL,
                           extra_facts = NULL, contrast_weights = NULL, p_shape_skip = NULL,
                           technical_panels = TRUE, gene_sets = NULL, p_floor = NULL) {
  n_features <- n_features %||% .default("n_features")
  fdr <- fdr %||% .default("fdr")
  effect_cutoff <- effect_cutoff %||% .default("effect_cutoff")
  p_floor <- p_floor %||% .default("p_floor")

  rows <- .rows(table)
  has_stat <- "stat" %in% names(rows)
  informative <- if (has_stat) rows[rows$stat != 0, , drop = FALSE] else rows
  ok <- .below(rows$p_adj, fdr) & abs(rows$effect) > min_effect & !.excluded(rows, exclude)
  eligible <- rows[ok, , drop = FALSE]
  up <- .rank(eligible[eligible$effect > 0, , drop = FALSE], has_stat)
  down <- .rank(eligible[eligible$effect < 0, , drop = FALSE], has_stat)
  top_up <- utils::head(up, n_features)
  top_down <- utils::head(down, n_features)

  h <- .p_histogram(informative$p)
  shape <- if (!is.null(p_shape_skip)) list(status = "not_run", note = p_shape_skip) else .p_shape(h)
  q <- if (nrow(informative)) stats::quantile(informative$effect, c(0.01, 0.25, 0.5, 0.75, 0.99), names = FALSE, type = 7) else rep(NA_real_, 5)
  facts <- list(
    fdr = fdr,
    effect_cutoff = effect_cutoff,
    n_hits = nrow(eligible),
    n_up = nrow(up),
    n_down = nrow(down),
    n_up_effect = sum(up$effect >= effect_cutoff),
    n_down_effect = sum(down$effect <= -effect_cutoff),
    n_near_miss = sum(!.below(rows$p_adj, fdr) & .below(rows$p_adj, 2 * fdr)),
    n_p_adj_missing = sum(is.na(rows$p_adj)),
    n_uninformative = nrow(rows) - nrow(informative),
    pi0 = .sig(.storey_pi0(informative$p), 3),
    p_histogram = h,
    effect_quantiles = list(q01 = .sig(q[1], 3), q25 = .sig(q[2], 3), q50 = .sig(q[3], 3),
                            q75 = .sig(q[4], 3), q99 = .sig(q[5], 3))
  )
  if (!isFALSE(technical_panels)) {
    is_sig <- .below(rows$p_adj, fdr) & abs(rows$effect) > min_effect
    significant <- .rank(rows[is_sig, , drop = FALSE], has_stat)
    tech <- .technical_signatures(.label(rows), .label(significant),
                                  if (isTRUE(technical_panels)) NULL else technical_panels)
    facts[names(tech)] <- tech
  }
  if (length(extra_facts)) facts[names(extra_facts)] <- extra_facts

  enrichment <- NULL
  if (!is.null(gene_sets)) {
    e_up <- .enrichment(.label(rows), .label(up), gene_sets)
    e_down <- .enrichment(.label(rows), .label(down), gene_sets)
    enrichment <- list(n_sets_tested = e_up$n_tested, up = e_up$sets, down = e_down$sets)
  }

  panel_state <- .empty_object()
  for (name in panel) {
    i <- match(name, rows$feature)
    if (is.na(i) && "symbol" %in% names(rows)) i <- match(name, rows$symbol)
    panel_state[[name]] <- if (is.na(i)) {
      NA  # not measured or filtered out: null, never 0
    } else {
      .drop_null(list(
        effect = .sig(rows$effect[i]),
        ave_expr = if ("ave_expr" %in% names(rows)) .sig(rows$ave_expr[i]),
        stat = if (has_stat) .sig(rows$stat[i]),
        p_adj = .sig(max(rows$p_adj[i], p_floor))
      ))
    }
  }

  rule <- sprintf("p_adj < %s", .fmt(fdr))
  if (min_effect > 0) rule <- paste0(rule, sprintf(" and |effect| > %s", .fmt(min_effect)))
  rule <- paste0(rule, sprintf(", ranked by %s, up to %d per direction",
                               if (has_stat) "|stat|" else "p then |effect|", as.integer(n_features)))

  .drop_null(list(
    unit_id = as.character(unit_id),
    unit_type = "contrast",
    definition = definition,
    contrast_weights = if (length(contrast_weights)) as.list(contrast_weights),
    size = c(list(n_tested = nrow(rows)), size),
    facts = facts,
    top_features = list(up = .feature_rows(top_up, p_floor), down = .feature_rows(top_down, p_floor)),
    panel = panel_state,
    tail = list(n_shown = nrow(top_up) + nrow(top_down), n_eligible = nrow(eligible), rule = rule),
    enrichment = enrichment,
    trust = list(.check("pvalue_distribution", shape$status, shape$note))
  ))
}

.cluster_unit <- function(table, unit_id, n_features = NULL, fdr = NULL, min_effect = 0, exclude = NULL,
                          panel_values = NULL, size = NULL, facts = NULL, technical_panels = TRUE,
                          gene_sets = NULL, p_floor = NULL) {
  n_features <- n_features %||% .default("n_features")
  fdr <- fdr %||% .default("fdr")
  p_floor <- p_floor %||% .default("p_floor")

  rows <- .rows(table)
  has_stat <- "stat" %in% names(rows)
  ok <- .below(rows$p_adj, fdr) & rows$effect > min_effect & !.excluded(rows, exclude)
  up <- .rank(rows[ok, , drop = FALSE], has_stat, signed = 1)
  top <- utils::head(up, n_features)
  rule <- sprintf("p_adj < %s and effect > %s, ranked by %s, top %d", .fmt(fdr), .fmt(min_effect),
                  if (has_stat) "stat" else "p then effect", as.integer(n_features))
  facts <- if (length(facts)) facts else .empty_object()
  if (!isFALSE(technical_panels)) {
    # scored before the name exclusion, so mitochondrial and ribosomal signatures stay visible
    is_sig <- .below(rows$p_adj, fdr) & rows$effect > min_effect
    significant <- .rank(rows[is_sig, , drop = FALSE], has_stat, signed = 1)
    tech <- .technical_signatures(.label(rows), .label(significant),
                                  if (isTRUE(technical_panels)) NULL else technical_panels)
    facts[names(tech)] <- tech
  }
  enrichment <- NULL
  if (!is.null(gene_sets)) {
    e_up <- .enrichment(.label(rows), .label(up), gene_sets)
    enrichment <- list(n_sets_tested = e_up$n_tested, up = e_up$sets)
  }
  .drop_null(list(
    unit_id = as.character(unit_id),
    unit_type = "cluster",
    size = if (length(size)) size else .empty_object(),
    facts = facts,
    top_features = .feature_rows(top, p_floor),
    panel = if (length(panel_values)) panel_values else .empty_object(),
    tail = list(n_shown = nrow(top), n_eligible = nrow(up), rule = rule),
    enrichment = enrichment
  ))
}

# Per-cluster facts from per-observation metadata, each compared with the whole dataset.
.cluster_facts <- function(labels, qc = NULL, doublet = NULL, batch = NULL, sample = NULL) {
  labels <- as.character(labels)
  n <- length(labels)
  clusters <- unique(labels)
  cfg <- llomics_contract()$batch_dominance

  reference <- list(definition = "medians and fractions across all observations in the dataset")
  for (key in names(qc)) reference[[paste0("median_", key)]] <- .sig(stats::median(qc[[key]], na.rm = TRUE))
  if (!is.null(doublet)) reference$pct_predicted_doublet <- .sig(100 * mean(as.logical(doublet)), 3)
  if (!is.null(batch)) {
    batch <- as.character(batch)
    batch_share <- table(batch) / n
  }

  facts <- list()
  flagged <- character(0)
  for (cl in clusters) {
    idx <- which(labels == cl)
    f <- list()
    for (key in names(qc)) {
      med <- stats::median(qc[[key]][idx], na.rm = TRUE)
      ref <- reference[[paste0("median_", key)]]
      f[[paste0("median_", key)]] <- .sig(med)
      f[[paste0(key, "_vs_reference")]] <- if (!is.na(ref) && ref != 0 && !is.na(med)) .sig(med / ref, 3) else NA
    }
    if (!is.null(doublet)) {
      pct <- 100 * mean(as.logical(doublet[idx]))
      ref <- reference$pct_predicted_doublet
      f$pct_predicted_doublet <- .sig(pct, 3)
      f$doublet_fraction_vs_reference <- if (ref != 0) .sig(pct / ref, 3) else NA
    }
    if (!is.null(sample)) f$n_samples_contributing <- length(unique(as.character(sample[idx])))
    if (!is.null(batch)) {
      names_sorted <- sort(unique(batch[idx]), method = "radix")
      counts <- vapply(names_sorted, function(b) sum(batch[idx] == b), numeric(1))
      top_batch <- names_sorted[which.max(counts)]
      share <- max(counts) / length(idx)
      f$dominant_batch <- top_batch
      f$dominant_batch_fraction <- .sig(share, 3)
      if (share >= cfg$cluster_fraction && batch_share[[top_batch]] < cfg$dataset_share_below) {
        flagged <- c(flagged, sprintf("%s (%.0f%% %s)", cl, 100 * share, top_batch))
      }
    }
    facts[[cl]] <- f
  }
  list(facts = facts, reference = reference, flagged = flagged)
}

.cluster_trust <- function(test, n_samples = NULL, doublet_note = NULL, n_batches = NULL, batch_flagged = character(0)) {
  cfg <- llomics_contract()$batch_dominance
  samples <- if (!is.null(n_samples)) sprintf("%d biological sample(s)", as.integer(n_samples)) else
    "the number of biological samples is unknown (no sample column supplied)"
  checks <- list(
    .check("pseudoreplication", "warn", sprintf(
      "%s with cells as replicates; %s. P-values rank markers and are not evidence of a population-level difference.",
      test, samples)),
    if (!is.null(doublet_note)) .check("doublet_detection", "pass", doublet_note) else
      .check("doublet_detection", "not_run", "no predicted-doublet column was supplied"),
    .check("cluster_stability", "not_run", "cannot be computed from one clustering; declare it with trust=")
  )
  batch_check <- if (is.null(n_batches)) {
    .check("batch_integration", "not_run", "no batch column supplied; declare a single batch with trust=")
  } else if (n_batches == 1) {
    .check("batch_integration", "not_applicable", "all observations come from one batch")
  } else if (length(batch_flagged)) {
    .check("batch_integration", "warn", paste0("clusters dominated by one batch: ", paste(batch_flagged, collapse = ", ")))
  } else {
    .check("batch_integration", "pass", sprintf(
      "no cluster is %.0f%% or more from a batch that holds under %.0f%% of observations",
      100 * cfg$cluster_fraction, 100 * cfg$dataset_share_below))
  }
  c(checks, list(batch_check))
}

# ---- payloads --------------------------------------------------------------------------------

.envelope <- function(method, run, dataset, trust, units) {
  contract <- llomics_contract()
  if (!method %in% names(contract$methods)) stop("unknown method '", method, "'")
  list(llomics = contract$version, method = method, run = run, dataset = dataset, trust = trust, units = units)
}

.standard_filters <- function(fdr, n_features, min_effect, exclude, ranked_by, extra = list()) {
  filters <- extra
  if (!is.null(exclude)) {
    filters <- c(filters, list(list(rule = paste0("top features exclude names matching ", exclude),
                                    why = "technical signal; applied to the feature list only, not to the test")))
  }
  filters <- c(filters, list(list(rule = sprintf("p_adj < %s", .fmt(fdr)), why = "FDR control")))
  if (min_effect > 0) {
    filters <- c(filters, list(list(rule = sprintf("|effect| > %s", .fmt(min_effect)), why = "minimum effect size")))
  }
  c(filters, list(list(rule = sprintf("top %d by %s; never padded", as.integer(n_features), ranked_by), why = "payload size")))
}

#' Build a differential-analysis payload from standard-column tables
#'
#' The generic entry point for any tool: one data.frame per contrast with the columns
#' `feature`, `effect`, `p`, `p_adj` and optionally `symbol`, `stat`, `ave_expr`, `pct_in`, `pct_out`.
#'
#' @param tables Named list of data.frames, one per contrast.
#' @param definitions Named character vector: each contrast in words, including the direction.
#' @param n_features How many features to load per direction. Never padded.
#' @param panel Features reported for every contrast whether or not they rank.
#' @param group_sizes Named vector of independent replicates per group.
#' @param batch_crosstab Named list by group of named lists by batch with observation counts.
#' @param trust Named list of declarations, e.g. `list(batch_confounding = c("not_applicable", "one run"))`.
#' @param technical_panels TRUE scores the shipped artefact panels against the hits, FALSE skips it,
#'   a named list replaces the panels.
#' @param gene_sets Named list of character vectors (see `llomics_read_gmt()`) for over-representation
#'   of the up and down hits. `gene_sets_name` records where they came from.
#' @export
llomics_dge_payload <- function(tables, definitions, analysis_id, question, tool, steps, n_observations,
                                observation_type = "sample", group_sizes = NULL, units_are_replicates = TRUE,
                                batch_crosstab = NULL, batch_in_design = NA, batch_adjustable = TRUE,
                                single_batch = FALSE, mean_variance = NULL, mean_variance_detected = NULL,
                                correction = "Benjamini-Hochberg", min_replicates = NULL, n_features = NULL,
                                fdr = NULL, min_effect = 0, effect_cutoff = NULL, exclude = NULL,
                                panel = character(0), sizes = NULL, contrast_weights = NULL, filters = list(),
                                source_tables = NULL, run_extra = NULL, dataset_reference = NULL,
                                p_shape_skip = NULL, technical_panels = TRUE, gene_sets = NULL,
                                gene_sets_name = NULL, trust = NULL, p_floor = NULL) {
  n_features <- n_features %||% .default("n_features")
  fdr <- fdr %||% .default("fdr")
  p_floor <- p_floor %||% .default("p_floor")
  missing <- setdiff(names(tables), names(definitions))
  if (length(missing)) {
    stop("`definitions` needs an entry in words, including the direction, for: ", paste(missing, collapse = ", "))
  }

  units <- lapply(names(tables), function(uid) {
    .contrast_unit(tables[[uid]], unit_id = uid, definition = definitions[[uid]], n_features = n_features,
                   fdr = fdr, min_effect = min_effect, effect_cutoff = effect_cutoff, exclude = exclude,
                   panel = panel, size = sizes[[uid]], contrast_weights = contrast_weights[[uid]],
                   p_shape_skip = p_shape_skip, technical_panels = technical_panels, gene_sets = gene_sets,
                   p_floor = p_floor)
  })
  shape_status <- vapply(units, function(u) u$trust[[1]]$status, character(1))
  shape_note <- vapply(units, function(u) sprintf("%s: %s", u$unit_id, u$trust[[1]]$note), character(1))
  has_stat <- all(vapply(tables, function(t) "stat" %in% names(t), logical(1)))

  checks <- list(
    .replicates_check(group_sizes, min_replicates, units_are_replicates),
    if (!is.null(p_shape_skip)) .check("pvalue_distribution", "not_run", p_shape_skip) else
      .check("pvalue_distribution", .worst(shape_status), paste(shape_note, collapse = " | ")),
    .batch_confounding_check(batch_crosstab, batch_in_design, batch_adjustable, single_batch),
    .multiple_testing_check(length(units), correction),
    .mean_variance_check(mean_variance, mean_variance_detected)
  )
  run <- c(
    list(
      analysis_id = analysis_id,
      question = question,
      tool = tool,
      steps = I(as.character(steps)),
      filters = .standard_filters(fdr, n_features, min_effect, exclude,
                                  if (has_stat) "|stat| per direction" else "p then |effect| per direction", filters),
      p_floor = p_floor,
      n_p_values_floored = sum(vapply(tables, .count_floored, numeric(1), p_floor = p_floor))
    ),
    .analysis_notes(technical_panels, gene_sets, gene_sets_name),
    run_extra,
    list(source_tables = if (length(source_tables)) as.list(source_tables) else .empty_object())
  )
  dataset <- list(
    n_observations = n_observations,
    observation_type = observation_type,
    n_features = max(vapply(tables, nrow, numeric(1))),
    reference = c(list(n_contrasts = length(units)), dataset_reference)
  )
  .envelope("dge", run, dataset, .apply_trust_overrides(checks, trust), units)
}

#' Build a clusters payload from standard-column marker tables and per-observation metadata
#'
#' @param tables Named list of data.frames, one per cluster (markers against the rest).
#' @param labels Cluster label of every observation.
#' @param qc Named list of numeric QC vectors, one value per observation.
#' @param doublet,batch,sample Per-observation predicted-doublet flag, batch and biological sample.
#' @param panel_values Named list by cluster of named lists of panel values (NA for not measured).
#' @export
llomics_clusters_payload <- function(tables, analysis_id, question, tool, steps, test, labels, qc = NULL,
                                     doublet = NULL, batch = NULL, sample = NULL, panel_values = NULL,
                                     extra_facts = NULL, n_dataset_features = NULL, observation_type = "cell",
                                     n_features = NULL, fdr = NULL, min_effect = 0, exclude = NULL,
                                     filters = list(), source_tables = NULL, run_extra = NULL,
                                     technical_panels = TRUE, gene_sets = NULL, gene_sets_name = NULL,
                                     trust = NULL, p_floor = NULL) {
  n_features <- n_features %||% .default("n_features")
  fdr <- fdr %||% .default("fdr")
  p_floor <- p_floor %||% .default("p_floor")
  labels <- as.character(labels)
  n <- length(labels)
  cf <- .cluster_facts(labels, qc = qc, doublet = doublet, batch = batch, sample = sample)

  units <- lapply(names(tables), function(uid) {
    f <- cf$facts[[uid]] %||% list()
    extra <- extra_facts[[uid]]
    if (length(extra)) f[names(extra)] <- extra
    n_cells <- sum(labels == uid)
    .cluster_unit(tables[[uid]], unit_id = uid, n_features = n_features, fdr = fdr, min_effect = min_effect,
                  exclude = exclude, panel_values = panel_values[[uid]], facts = f,
                  technical_panels = technical_panels, gene_sets = gene_sets, p_floor = p_floor,
                  size = list(n_cells = n_cells, pct_of_total = .sig(100 * n_cells / n, 3)))
  })

  doublet_note <- if (!is.null(doublet)) {
    sprintf("%s%% of observations predicted doublets", .fmt(cf$reference$pct_predicted_doublet))
  }
  checks <- .cluster_trust(
    test = test,
    n_samples = if (!is.null(sample)) length(unique(as.character(sample))),
    doublet_note = doublet_note,
    n_batches = if (!is.null(batch)) length(unique(as.character(batch))),
    batch_flagged = cf$flagged)
  has_stat <- all(vapply(tables, function(t) "stat" %in% names(t), logical(1)))
  run <- c(
    list(
      analysis_id = analysis_id,
      question = question,
      tool = tool,
      steps = I(as.character(steps)),
      filters = .standard_filters(fdr, n_features, min_effect, exclude, if (has_stat) "stat" else "p then effect", filters),
      p_floor = p_floor,
      n_p_values_floored = sum(vapply(tables, .count_floored, numeric(1), p_floor = p_floor)),
      panel_units = "mean expression across the unit's observations, in the units of the tested matrix"
    ),
    .analysis_notes(technical_panels, gene_sets, gene_sets_name),
    run_extra,
    list(source_tables = if (length(source_tables)) as.list(source_tables) else .empty_object())
  )
  dataset <- list(
    n_observations = n,
    observation_type = observation_type,
    n_biological_samples = if (!is.null(sample)) length(unique(as.character(sample))) else NA,
    n_features = n_dataset_features %||% max(vapply(tables, nrow, numeric(1))),
    reference = cf$reference
  )
  .envelope("clusters", run, dataset, .apply_trust_overrides(checks, trust), units)
}

#' Write a payload to a JSON file
#' @export
llomics_write <- function(payload, out) {
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, out, auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE)
  message("llomics: wrote ", out, " (", length(payload$units), " unit", if (length(payload$units) != 1) "s", ")")
  invisible(out)
}

#' Validate a payload file against the contract
#'
#' Runs the shared Python validator shipped with the package. Needs a `python` (3.9 or later) on the path.
#' @export
llomics_validate <- function(path, python = Sys.which(c("python3", "python"))) {
  python <- python[nzchar(python)][1]
  if (is.na(python)) stop("no python found on the path; run validate_payload.py by hand")
  script <- system.file("validate_payload.py", package = "llomics")
  status <- system2(python, c(shQuote(script), shQuote(path)))
  invisible(status == 0)
}
