# Regenerate the benchmark payloads. Needs limma and the llomics R package. From the repository root:
#
#   Rscript benchmark/make_scenarios.R
#
# Each scenario simulates log-expression for 12 samples (unless stated) over real human gene
# symbols, fits limma, and writes benchmark/payloads/<id>.llomics.json through llomics_limma().
# The truth of each scenario and the expected decision are in benchmark/scenarios.json.

suppressPackageStartupMessages({library(limma); library(llomics)})
set.seed(20261117)

ISG <- c("IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "ISG20", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL",
         "IFI6", "IFI27", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFITM1", "IFITM3", "RSAD2", "CXCL10", "STAT1",
         "STAT2", "IRF7", "IRF9", "HERC5", "HERC6", "USP18", "DDX58", "DDX60", "EPSTI1", "XAF1", "SAMD9",
         "SAMD9L", "CMPK2", "GBP1", "GBP4", "BST2", "LY6E", "PARP9", "PARP14", "TRIM22", "SP110", "EIF2AK2",
         "PLSCR1", "LGALS3BP", "SIGLEC1", "HELZ2", "ZBP1", "DTX3L", "NMI", "TRIM21", "ADAR", "UBE2L6")
HYPOXIA <- c("VEGFA", "CA9", "ADM", "BNIP3", "BNIP3L", "EGLN3", "EGLN1", "NDRG1", "SLC2A1", "SLC2A3", "PGK1",
             "LDHA", "ANKRD37", "PDK1", "P4HA1", "P4HA2", "ENO1", "ALDOA", "HK2", "PFKFB3", "PFKFB4", "KDM3A",
             "LOX", "BHLHE40", "DDIT4", "ANGPTL4", "ERO1A", "MXI1", "STC2", "HILPDA", "PFKP", "TPI1", "GAPDH",
             "PGAM1", "INSIG2", "FAM162A", "CA12", "PLOD2", "AK4", "KDM5B")
SEX_F <- c("XIST", "TSIX", "JPX")
SEX_M <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY", "TXLNGY", "NLGN4Y", "TMSB4Y", "TTTY14",
           "TTTY15", "PRKY", "TBL1Y")
SEX_ESCAPE <- c("KDM5C", "KDM6A", "DDX3X", "EIF1AX", "ZFX", "PUDP", "SMC1A", "CA5B", "EIF2S3", "TXLNG")
STRESS <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "IER5", "ATF3",
            "DUSP1", "DUSP2", "NR4A1", "NR4A2", "NR4A3", "KLF2", "KLF4", "KLF6", "HSPA1A", "HSPA1B", "HSPA8",
            "HSPH1", "DNAJB1", "DNAJA1", "HSP90AA1", "HSPB1", "HSPE1", "ZFP36", "BTG1", "BTG2", "SOCS3", "MYC",
            "CCN1", "PPP1R15A", "MCL1", "NFKBIA", "CEBPB", "CEBPD", "GADD45B", "SGK1", "TOB1", "UBC")
LYMPHOCYTE <- c("PTPRC", "CD3E", "CD3D", "CD3G", "CD2", "CD5", "CD7", "CD19", "MS4A1", "CD79A", "CD79B", "CD4",
                "CD8A", "CD8B", "LCK", "ZAP70", "CCL5", "GZMK", "NKG7", "CD52", "CORO1A", "LAPTM5", "CD37",
                "CD53", "IKZF1", "ITGB2", "CD48", "RAC2", "ARHGDIB", "LCP1", "CXCR4", "CCR7", "SELL", "IL7R",
                "TRAC", "CD247", "CD69", "LTB", "TRBC1", "TRBC2", "CD27", "CD28", "ITK", "IL2RG")

background <- readLines("benchmark/data/gene_symbols.txt")
genes <- unique(c(background, ISG, HYPOXIA, SEX_F, SEX_M, SEX_ESCAPE, STRESS, LYMPHOCYTE))
n_genes <- length(genes)
baseline <- stats::setNames(runif(n_genes, 3, 11), genes)
baseline[c(SEX_F, SEX_M)] <- 1                      # off unless the sex expresses them
gene_sd <- stats::setNames(sqrt(0.09 * 8 / rchisq(n_genes, df = 8)), genes)

`%||%` <- function(a, b) if (is.null(a)) b else a
simulate <- function(targets, bio = NULL) {
  n <- nrow(targets)
  bio <- bio %||% seq_len(n)                        # which biological sample each column measures
  bio_values <- matrix(rnorm(n_genes * max(bio), mean = baseline, sd = gene_sd), nrow = n_genes)
  expr <- bio_values[, bio, drop = FALSE]
  if (anyDuplicated(bio)) expr <- expr + matrix(rnorm(n_genes * n, sd = 0.08), nrow = n_genes)
  dimnames(expr) <- list(genes, targets$sample)
  expr
}
shift <- function(expr, set, cols, low = 1.5, high = 3) {
  set <- intersect(set, rownames(expr))
  expr[set, cols] <- expr[set, cols] + runif(length(set), low, high)
  expr
}
make_targets <- function(n_per_group = 6, aliased = FALSE) {
  group <- rep(c("control", "treated"), each = n_per_group)
  batch <- if (aliased) rep(c("b1", "b2"), each = n_per_group) else rep(c("b1", "b2"), length.out = 2 * n_per_group)
  data.frame(sample = sprintf("S%02d", seq_along(group)), group = group, batch = batch)
}

run <- function(id, expr, targets, definition, question, panel, batch_in_design = TRUE, give_batch = TRUE) {
  group <- factor(targets$group)
  design <- if (batch_in_design) model.matrix(~ 0 + group + targets$batch) else model.matrix(~ 0 + group)
  colnames(design)[1:2] <- c("control", "treated")
  if (ncol(design) > 2) colnames(design)[3] <- "batch_b2"
  fit <- lmFit(expr, design)
  fit <- contrasts.fit(fit, makeContrasts(treated_vs_control = treated - control, levels = design))
  fit <- eBayes(fit, trend = TRUE)
  payload <- llomics_limma(
    fit, coefs = "treated_vs_control", definitions = c(treated_vs_control = definition),
    targets = targets, group = "group", batch = if (give_batch) "batch", batch_in_design = batch_in_design,
    analysis_id = id, question = question,
    steps = c("simulated log2 expression over human gene symbols",
              if (batch_in_design) "lmFit ~ 0 + group + batch" else "lmFit ~ 0 + group",
              "contrasts.fit treated - control", "eBayes(trend = TRUE)"),
    panel = panel, out_dir = NULL)
  path <- file.path("benchmark/payloads", paste0(id, ".llomics.json"))
  payload$run$source_tables <- list(note = "benchmark payload; the full table is not kept")
  suppressMessages(llomics_write(payload, path))
  u <- payload$units[[1]]
  cat(sprintf("%-22s hits %4d (up %3d, down %3d)\n", id, u$facts$n_hits, u$facts$n_up, u$facts$n_down))
}

HOUSEKEEPING <- c("ACTB", "GAPDH")
treated <- function(t) t$group == "treated"

# 1. clean_interferon
t <- make_targets(); e <- shift(simulate(t), ISG, treated(t))
run("clean_interferon", e, t, "interferon-beta 6 h - untreated, primary fibroblasts, adjusted for batch; positive effect = higher after interferon",
    "Which genes respond to interferon-beta?", c("ISG15", "MX1", HOUSEKEEPING))

# 2. clean_hypoxia
t <- make_targets(); e <- shift(simulate(t), HYPOXIA, treated(t), 1.2, 2.5)
run("clean_hypoxia", e, t, "1% oxygen 24 h - 21% oxygen, cultured cells, adjusted for batch; positive effect = higher in hypoxia",
    "Which genes respond to low oxygen?", c("VEGFA", "CA9", "ACTB"))

# 3. null
t <- make_targets(); e <- simulate(t)
run("null", e, t, "compound X 6 h - vehicle, adjusted for batch; positive effect = higher with compound X",
    "Which genes respond to compound X?", HOUSEKEEPING)

# 4. batch_aliased: a real response, plus a batch effect that cannot be separated from it
t <- make_targets(aliased = TRUE); e <- shift(simulate(t), ISG, treated(t))
e <- shift(e, sample(background, 600), t$batch == "b2", 0.8, 1.6)
run("batch_aliased", e, t, "interferon-beta 6 h - untreated, primary fibroblasts; positive effect = higher after interferon",
    "Which genes respond to interferon-beta?", c("ISG15", "MX1", HOUSEKEEPING), batch_in_design = FALSE)

# 5. no_batch_info: sound analysis, but the extractor is not told about batch
t <- make_targets(); e <- shift(simulate(t), ISG, treated(t))
run("no_batch_info", e, t, "interferon-beta 6 h - untreated, primary fibroblasts; positive effect = higher after interferon",
    "Which genes respond to interferon-beta?", c("ISG15", "MX1", HOUSEKEEPING), batch_in_design = TRUE, give_batch = FALSE)

# 5b. unmodelled_batch: balanced batches, strong batch effect, batch left out of the design
t <- make_targets(); e <- shift(simulate(t), ISG, treated(t))
e <- shift(e, sample(background, 1500), t$batch == "b2", 1.5, 2.5)
run("unmodelled_batch", e, t, "interferon-beta 6 h - untreated, primary fibroblasts; positive effect = higher after interferon",
    "Which genes respond to interferon-beta?", c("ISG15", "MX1", HOUSEKEEPING), batch_in_design = FALSE)

# 6. sex_imbalance: treated all female, control all male, sex not recorded
t <- make_targets(); e <- simulate(t)
e <- shift(e, SEX_F, treated(t), 6, 9); e <- shift(e, SEX_M, !treated(t), 5, 8)
e <- shift(e, SEX_ESCAPE, treated(t), 0.5, 0.9)
run("sex_imbalance", e, t, "drug 4 weeks - placebo, whole blood, adjusted for batch; positive effect = higher on drug",
    "Which genes respond to the drug?", HOUSEKEEPING)

# 7. stress_handling: treated samples processed later
t <- make_targets(); e <- shift(simulate(t), STRESS, treated(t), 1.2, 3)
run("stress_handling", e, t, "knockout - wild type, sorted liver macrophages, adjusted for batch; positive effect = higher in knockout",
    "Which genes depend on the knocked-out gene?", HOUSEKEEPING)

# 8. composition: more lymphocytes in treated tissue
t <- make_targets(); e <- shift(simulate(t), LYMPHOCYTE, treated(t), 1.2, 2.5)
run("composition", e, t, "high-fat diet 12 weeks - chow, bulk liver tissue, adjusted for batch; positive effect = higher on high-fat diet",
    "Which hepatocyte genes respond to the diet?", c("ALB", "APOA1", "ACTB"))

# 9. two_replicates: one batch, so the design is group only and hits do appear
t <- make_targets(n_per_group = 2); t$batch <- "b1"; e <- shift(simulate(t), ISG, treated(t))
run("two_replicates", e, t, "interferon-beta 6 h - untreated, primary fibroblasts; positive effect = higher after interferon",
    "Which genes respond to interferon-beta?", c("ISG15", "MX1", HOUSEKEEPING), batch_in_design = FALSE)

# 10. technical_replicates: 3 biological samples per group, each measured twice, treated as 6
t <- make_targets(); bio <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)
t$batch <- rep(c("b1", "b2"), times = 6)
e <- simulate(t, bio = bio)
run("technical_replicates", e, t, "compound Y 6 h - vehicle, adjusted for batch; positive effect = higher with compound Y",
    "Which genes respond to compound Y?", HOUSEKEEPING)
