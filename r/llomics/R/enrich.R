# Technical-panel signatures and gene-set enrichment, computed in code.
# Mirrors python/llomics/enrich.py.

.panels_env <- new.env(parent = emptyenv())

#' The technical panels shipped with llomics (sex, stress, cell cycle, haemoglobin, mitochondrial, ribosomal)
#' @export
llomics_panels <- function() {
  if (is.null(.panels_env$panels)) {
    path <- system.file("llomics-panels.json", package = "llomics")
    if (!nzchar(path)) stop("llomics-panels.json not found: run python tools/sync.py, then reinstall the package")
    .panels_env$panels <- jsonlite::fromJSON(path, simplifyVector = FALSE)$panels
  }
  .panels_env$panels
}

#' Read gene sets from a GMT file
#'
#' One set per line: name, description, then members, tab-separated.
#' @return A named list of character vectors.
#' @export
llomics_read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  sets <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) >= 3 && nzchar(parts[1])) sets[[parts[1]]] <- parts[-(1:2)][nzchar(parts[-(1:2)])]
  }
  sets
}

# P(X >= k) when drawing n from N of which K are marked
.hypergeom_sf <- function(k, N, K, n) {
  if (k <= 0) return(1)
  if (k > min(K, n)) return(0)
  min(1, stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE))
}

.in_panel <- function(labels_upper, panel) {
  genes <- toupper(unlist(panel$genes))
  hit <- labels_upper %in% genes
  if (!is.null(panel$pattern)) hit <- hit | grepl(panel$pattern, labels_upper, ignore.case = TRUE, perl = TRUE)
  hit
}

# Facts scoring each technical panel against the hits. When no panel feature was tested at all,
# the flag count is NA: the question could not be asked, which is not a zero.
.technical_signatures <- function(tested_labels, ranked_hit_labels, panels = NULL) {
  panels <- panels %||% llomics_panels()
  cfg <- llomics_contract()$technical_flag
  tested <- toupper(as.character(tested_labels))
  hits <- toupper(as.character(ranked_hit_labels))
  top <- utils::head(hits, cfg$top_n)
  signatures <- .empty_object()
  flags <- character(0)
  measured_total <- 0L
  for (name in names(panels)) {
    panel <- panels[[name]]
    n_measured <- sum(.in_panel(tested, panel))
    measured_total <- measured_total + n_measured
    n_hits <- sum(.in_panel(hits, panel))
    if (n_hits == 0) next
    top_share <- sum(.in_panel(top, panel)) / length(top)
    p <- max(.hypergeom_sf(n_hits, length(tested), n_measured, length(hits)), .default("p_floor"))
    signatures[[name]] <- list(n_measured = n_measured, n_hits = n_hits, top_share = signif(top_share, 3), p = signif(p, 4))
    if (n_hits >= cfg$min_hits && p < cfg$max_p && top_share >= cfg$min_top_share) flags <- c(flags, name)
  }
  list(
    n_panel_features_measured = measured_total,
    technical_signatures = signatures,
    technical_flags = I(flags),
    n_technical_flags = if (measured_total > 0) length(flags) else NA
  )
}

# Over-representation of gene sets among the hits, against the tested features. Returns the
# significant sets (best first, capped) and the number of sets tested.
.enrichment <- function(tested_labels, ranked_hit_labels, gene_sets) {
  cfg <- llomics_contract()$enrichment
  universe <- unique(toupper(as.character(tested_labels)))
  hit_original <- as.character(ranked_hit_labels)
  keep <- !duplicated(toupper(hit_original))
  hit_original <- hit_original[keep]
  hit_upper <- toupper(hit_original)

  tested <- list()
  for (name in names(gene_sets)) {
    inside <- intersect(unique(toupper(as.character(unlist(gene_sets[[name]])))), universe)
    if (length(inside) < cfg$min_set_size) next
    overlap <- hit_original[hit_upper %in% inside]
    tested[[length(tested) + 1L]] <- list(
      name = name, n_set = length(inside), overlap = overlap,
      p = .hypergeom_sf(length(overlap), length(universe), length(inside), length(hit_upper)))
  }
  if (!length(tested)) return(list(sets = list(), n_tested = 0L))

  p <- vapply(tested, function(t) t$p, numeric(1))
  p_adj <- stats::p.adjust(p, method = "BH")
  names_all <- vapply(tested, function(t) t$name, character(1))
  n_overlap <- vapply(tested, function(t) length(t$overlap), numeric(1))
  keep <- which(p_adj < cfg$fdr & n_overlap >= cfg$min_overlap)
  keep <- keep[order(p[keep], names_all[keep], method = "radix")]
  keep <- utils::head(keep, cfg$max_sets)
  sets <- lapply(keep, function(i) {
    t <- tested[[i]]
    list(set = t$name, n_set = t$n_set, n_overlap = length(t$overlap),
         expected = signif(length(hit_upper) * t$n_set / length(universe), 3),
         p_adj = signif(max(p_adj[i], .default("p_floor")), 4),
         genes = I(utils::head(t$overlap, cfg$max_genes)))
  })
  list(sets = sets, n_tested = length(tested))
}

.analysis_notes <- function(technical_panels, gene_sets, gene_sets_name) {
  contract <- llomics_contract()
  cfg_t <- contract$technical_flag
  cfg_e <- contract$enrichment
  notes <- list()
  if (!isFALSE(technical_panels)) {
    notes$technical_panels <- sprintf(
      "each panel is scored against all significant features; a panel is flagged at %d or more hits, hypergeometric p < %s, and at least %.0f%% of the top %d hits",
      as.integer(cfg_t$min_hits), .fmt(cfg_t$max_p), 100 * cfg_t$min_top_share, as.integer(cfg_t$top_n))
  }
  if (!is.null(gene_sets)) {
    notes$enrichment <- list(
      gene_sets = gene_sets_name %||% "supplied by the analyst",
      n_sets = length(gene_sets),
      method = sprintf(
        "one-sided hypergeometric test of the hits against the tested features, sets with at least %d tested members, Benjamini-Hochberg across sets, p_adj < %s",
        as.integer(cfg_e$min_set_size), .fmt(cfg_e$fdr)))
  }
  notes
}
