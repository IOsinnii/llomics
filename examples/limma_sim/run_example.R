# Smoke test for the limma adapter on simulated data. Install the R package first
# (R CMD INSTALL r/llomics), then run from the repository root:
#
#   Rscript examples/limma_sim/run_example.R [out_dir]
#
# Simulates 5,000 features x 12 samples (6 control, 6 stimulated, 2 balanced batches) with
# 200 features up and 120 down, fits limma, and writes an llomics payload to out_dir
# (default: examples/limma_sim/out).

suppressPackageStartupMessages(library(limma))
library(llomics)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else "examples/limma_sim/out"

set.seed(1)
n_genes <- 5000
targets <- data.frame(
  sample = sprintf("S%02d", 1:12),
  group = rep(c("control", "stimulated"), each = 6),
  batch = rep(c("b1", "b2"), times = 6)
)

sd_gene <- sqrt(0.25 * 8 / rchisq(n_genes, df = 8))
expr <- matrix(rnorm(n_genes * 12, mean = rep(runif(n_genes, 3, 11), 12), sd = rep(sd_gene, 12)), nrow = n_genes)
rownames(expr) <- sprintf("G%04d", seq_len(n_genes))
colnames(expr) <- targets$sample

stim <- targets$group == "stimulated"
expr[1:200, stim] <- expr[1:200, stim] + 1.5
expr[201:320, stim] <- expr[201:320, stim] - 1.2
batch_genes <- sample(n_genes, 500)
expr[batch_genes, targets$batch == "b2"] <- expr[batch_genes, targets$batch == "b2"] + 0.8

group <- factor(targets$group)
batch <- factor(targets$batch)
design <- model.matrix(~ 0 + group + batch)
colnames(design) <- c("control", "stimulated", "batch_b2")

fit <- lmFit(expr, design)
fit <- contrasts.fit(fit, makeContrasts(stim_vs_ctrl = stimulated - control, levels = design))
fit <- eBayes(fit, trend = TRUE)
fit$genes <- data.frame(symbol = sprintf("SYM%04d", seq_len(n_genes)))

llomics_limma(
  fit,
  coefs = "stim_vs_ctrl",
  definitions = c(stim_vs_ctrl = "stimulated - control, adjusted for batch; positive effect = higher in stimulated"),
  targets = targets, group = "group", batch = "batch", batch_in_design = TRUE,
  analysis_id = "limma-sim",
  question = "Which features respond to stimulation? (simulated data: features have no biological meaning)",
  steps = c("simulated log-expression, 5000 features x 12 samples",
            "lmFit ~ 0 + group + batch", "contrasts.fit stimulated - control", "eBayes(trend = TRUE)"),
  panel = c("G0001", "SYM4000", "NOT_MEASURED"),
  symbol_col = "symbol",
  out_dir = out_dir
)
