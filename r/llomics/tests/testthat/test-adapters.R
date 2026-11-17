# The core is compared with the Python core by tests/test_parity.py at the repository root.
# These tests cover what only exists in R: the limma and Seurat adapters, and JSON writing.

check_status <- function(payload, name) {
  for (chk in payload$trust) if (chk$check == name) return(chk$status)
  NA_character_
}

written <- function(payload) {
  path <- tempfile(fileext = ".json")
  suppressMessages(llomics_write(payload, path))
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

sim_limma <- function(aliased = FALSE) {
  set.seed(1)
  n <- 800
  targets <- data.frame(group = rep(c("control", "stimulated"), each = 4),
                        batch = if (aliased) rep(c("b1", "b2"), each = 4) else rep(c("b1", "b2"), times = 4))
  expr <- matrix(rnorm(n * 8, mean = 7, sd = 0.4), nrow = n, dimnames = list(sprintf("G%03d", 1:n), NULL))
  expr[1:40, targets$group == "stimulated"] <- expr[1:40, targets$group == "stimulated"] + 2
  design <- stats::model.matrix(~ 0 + group, data = targets)
  colnames(design) <- c("control", "stimulated")
  fit <- limma::lmFit(expr, design)
  fit <- limma::contrasts.fit(fit, limma::makeContrasts(stim_vs_ctrl = stimulated - control, levels = design))
  fit <- limma::eBayes(fit, trend = TRUE)
  fit$genes <- data.frame(symbol = sprintf("SYM%03d", 1:n))
  list(fit = fit, targets = targets)
}

test_that("limma adapter builds a complete dge payload", {
  skip_if_not_installed("limma")
  s <- sim_limma()
  payload <- llomics_limma(
    s$fit, coefs = "stim_vs_ctrl",
    definitions = c(stim_vs_ctrl = "stimulated - control; positive effect = higher in stimulated"),
    targets = s$targets, group = "group", batch = "batch", batch_in_design = FALSE,
    analysis_id = "t", question = "q", steps = "lmFit", n_features = 5,
    panel = c("G001", "SYM700", "NOT_MEASURED"), symbol_col = "symbol", out_dir = NULL)
  expect_equal(payload$method, "dge")
  unit <- payload$units[[1]]
  expect_equal(nrow(unit$top_features$up), 5)
  expect_true(all(unit$top_features$up$feature %in% sprintf("G%03d", 1:40)))
  expect_equal(unit$contrast_weights, list(control = -1, stimulated = 1))
  expect_true(is.na(unit$panel$NOT_MEASURED))
  expect_false(is.null(unit$panel$SYM700$p_adj))
  expect_equal(check_status(payload, "replicates"), "pass")
  expect_equal(check_status(payload, "mean_variance"), "pass")      # eBayes(trend = TRUE) detected
  expect_equal(check_status(payload, "batch_confounding"), "warn")  # batch not in the design
  expect_equal(written(payload)$units[[1]]$panel$NOT_MEASURED, NULL) # written as JSON null
})

test_that("limma adapter fails the batch check when group and batch are aliased", {
  skip_if_not_installed("limma")
  s <- sim_limma(aliased = TRUE)
  payload <- llomics_limma(
    s$fit, coefs = "stim_vs_ctrl", definitions = c(stim_vs_ctrl = "stimulated - control"),
    targets = s$targets, group = "group", batch = "batch", batch_in_design = TRUE,
    analysis_id = "t", question = "q", steps = "lmFit", out_dir = NULL)
  expect_equal(check_status(payload, "batch_confounding"), "fail")
})

test_that("limma adapter refuses a fit without moderated statistics or a definition", {
  skip_if_not_installed("limma")
  s <- sim_limma()
  expect_error(llomics_limma(s$fit, coefs = "stim_vs_ctrl", definitions = c(other = "x"), targets = s$targets,
                             group = "group", analysis_id = "t", question = "q", steps = "s", out_dir = NULL),
               "definitions")
  raw <- s$fit
  raw$t <- NULL
  expect_error(llomics_limma(raw, coefs = "stim_vs_ctrl", definitions = c(stim_vs_ctrl = "x"), targets = s$targets,
                             group = "group", analysis_id = "t", question = "q", steps = "s", out_dir = NULL),
               "eBayes")
})

test_that("JSON keeps arrays as arrays and empty objects as objects", {
  table <- data.frame(feature = c("a", "b"), effect = c(2, -0.1), p = c(1e-6, 0.6), p_adj = c(1e-5, 0.6))
  payload <- llomics_dge_payload(list(x = table), c(x = "x - y"), analysis_id = "t", question = "q",
                                 tool = "t", steps = "one step", n_observations = 4)
  json <- written(payload)
  expect_type(json$run$steps, "list")
  expect_length(json$run$steps, 1)
  expect_length(json$units[[1]]$top_features$up, 1)       # one row is still an array of rows
  expect_length(json$units[[1]]$top_features$down, 0)
  expect_equal(json$units[[1]]$panel, structure(list(), names = character(0)))
  expect_equal(check_status(payload, "replicates"), "not_run")
})

test_that("Seurat FindAllMarkers wrapper builds a clusters payload", {
  skip_if_not_installed("Seurat")
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$donor <- rep(c("d1", "d2"), length.out = ncol(pbmc_small))
  payload <- llomics_FindAllMarkers(pbmc_small, n_features = 6, panel = c("CD3D", "MS4A1", "NOT_A_GENE"),
                                    question = "q", sample_col = "donor", batch_col = "groups", verbose = FALSE)
  expect_equal(payload$method, "clusters")
  expect_equal(vapply(payload$units, function(u) u$unit_id, ""), levels(Seurat::Idents(pbmc_small)))
  unit <- payload$units[[1]]
  expect_lte(nrow(unit$top_features), 6)
  expect_equal(unit$tail$n_shown, nrow(unit$top_features))
  expect_true(all(c("pct_in", "pct_out") %in% names(unit$top_features)))
  expect_true(is.na(unit$panel$NOT_A_GENE))
  expect_true(is.numeric(unit$panel$CD3D))
  expect_false(is.null(unit$facts$nearest_cluster))
  expect_equal(check_status(payload, "pseudoreplication"), "warn")
  expect_equal(check_status(payload, "doublet_detection"), "not_run")
  expect_equal(payload$dataset$n_biological_samples, 2)
  expect_error(llomics_FindAllMarkers(pbmc_small, batch_col = "nope", markers = data.frame(cluster = "0", gene = "a",
               avg_log2FC = 1, p_val = 0.01, p_val_adj = 0.01, pct.1 = 1, pct.2 = 0)), "not a metadata column")
})

test_that("Seurat FindMarkers wrapper marks a pre-filtered histogram as not judged", {
  skip_if_not_installed("Seurat")
  data("pbmc_small", package = "SeuratObject", envir = environment())
  pbmc_small$donor <- rep(c("d1", "d2"), length.out = ncol(pbmc_small))
  filtered <- llomics_FindMarkers(pbmc_small, ident.1 = "0", ident.2 = "1", n_features = 4,
                                  panel = c("LYZ", "NOT_A_GENE"), sample_col = "donor", verbose = FALSE)
  unit <- filtered$units[[1]]
  expect_equal(unit$unit_id, "0_vs_1")
  expect_match(unit$definition, "higher in 0")
  expect_true(all(unit$top_features$up$effect > 0))
  expect_true(all(unit$top_features$down$effect < 0))
  expect_true(is.na(unit$panel$NOT_A_GENE))
  expect_equal(check_status(filtered, "pvalue_distribution"), "not_run")
  expect_equal(check_status(filtered, "replicates"), "warn")
  expect_equal(check_status(filtered, "mean_variance"), "not_applicable")

  full <- llomics_FindMarkers(pbmc_small, ident.1 = "0", ident.2 = "1", logfc.threshold = 0, min.pct = 0, verbose = FALSE)
  expect_true(check_status(full, "pvalue_distribution") %in% c("pass", "warn"))
  expect_equal(full$units[[1]]$size$n_tested, nrow(pbmc_small))
})
