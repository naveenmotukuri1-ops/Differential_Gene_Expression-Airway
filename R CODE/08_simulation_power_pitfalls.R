# =============================================================================
# MODULE 08  --  Know your tools: null data, simulation, power, and the
#                classic mistakes (airway-based)
# =============================================================================
# Learning goals
#   * A "mock" comparison where nothing should be DE: does your pipeline agree?
#   * Simulate counts with KNOWN truth from real means/dispersions
#   * Measure real FDR and sensitivity of DESeq2 / edgeR / limma-voom vs n
#   * How many replicates do you need? (power curve)
#   * Pitfall 1: technical replicates treated as biological (pseudoreplication)
#   * Pitfall 2: double dipping -- heatmaps of DE genes "separate" pure noise
#
# Run:  Rscript 08_simulation_power_pitfalls.R   Outputs: results/08_simulation/
# Runtime: a few minutes. Needs results from 04.
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(pheatmap)
})
od <- out_dir("08_simulation")
set.seed(2024)

aw <- load_airway()
de2 <- readRDS(file.path(RESULTS_DIR, "04_deseq2", "airway_deseq2.rds"))
dds_fit <- de2$dds

# -----------------------------------------------------------------------------
banner("1. Mock comparison: control vs control")
# -----------------------------------------------------------------------------
ctrl <- aw$coldata$dex == "control"
cts_c <- aw$counts[, ctrl]; cd_c <- droplevels(aw$coldata[ctrl, ])
cts_c <- cts_c[rowSums(cts_c >= 10) >= 2, ]
cd_c$fake <- factor(c("A", "A", "B", "B"))     # arbitrary split of 4 donors
print(cd_c[, c("cell", "fake")])
dds_m <- DESeq(DESeqDataSetFromMatrix(cts_c, cd_c, ~ fake), quiet = TRUE)
res_m <- results(dds_m)
cat("Genes 'DE' at FDR 5% in a comparison with no real biology:",
    sum(res_m$padj < .05, na.rm = TRUE), "\n")
save_plot(ggplot(as.data.frame(res_m[res_m$baseMean > 10, ]), aes(pvalue)) +
            geom_histogram(breaks = seq(0, 1, .025)) +
            labs(title = "Mock comparison p-values", subtitle = "Should be ~uniform"),
          file.path(od, "mock_pvalues.png"))
cat("A handful of hits here are donor-specific genes (any 2-vs-2 split of donors\n",
    "differs somewhere). Run a mock comparison on your own data whenever you can.\n")

# -----------------------------------------------------------------------------
banner("2. A simulator built from real data")
# -----------------------------------------------------------------------------
# Take per-gene mean (control level) and dispersion from the airway fit.
mu0  <- rowMeans(counts(dds_fit, normalized = TRUE)[, dds_fit$dex == "control"])
disp <- dispersions(dds_fit)
ok   <- is.finite(disp) & mu0 > 5
pool <- data.frame(mu = mu0[ok], disp = disp[ok])
cat("Gene pool for simulation:", nrow(pool), "genes\n")

# disp_mult scales biological variability: 1 = this cell-line experiment
# (BCV ~0.1); 9 gives BCV ~0.3, typical of animal tissue or human samples.
simulate <- function(n_per_group, n_genes = 4000, frac_de = 0.1, sd_lfc = 1, disp_mult = 1) {
  g <- pool[sample(nrow(pool), n_genes, replace = TRUE), ]
  true_lfc <- numeric(n_genes)
  de <- sample(n_genes, round(frac_de * n_genes))
  true_lfc[de] <- sign(rnorm(length(de))) * (0.5 + abs(rnorm(length(de), 0, sd_lfc)))
  grp <- factor(rep(c("ctrl", "trt"), each = n_per_group))
  sf  <- exp(rnorm(2 * n_per_group, 0, 0.2))                 # unequal depth
  mu  <- outer(g$mu, sf) * 2^(outer(true_lfc, as.numeric(grp == "trt")))
  cts <- matrix(rnbinom(length(mu), mu = mu, size = 1 / (g$disp * disp_mult)), nrow = n_genes)
  storage.mode(cts) <- "integer"
  rownames(cts) <- paste0("g", seq_len(n_genes)); colnames(cts) <- paste0("s", seq_along(grp))
  list(counts = cts, group = grp, truth = true_lfc != 0, true_lfc = true_lfc)
}

run_methods <- function(sim) {
  cd <- data.frame(group = sim$group, row.names = colnames(sim$counts))
  keep <- rowSums(sim$counts >= 5) >= 2
  out <- list()
  # DESeq2
  d <- DESeq(DESeqDataSetFromMatrix(sim$counts[keep, ], cd, ~ group), quiet = TRUE)
  p <- rep(NA, length(keep)); p[keep] <- results(d)$padj; out$DESeq2 <- p
  # edgeR QL
  design <- model.matrix(~ group, cd)
  y <- normLibSizes(DGEList(sim$counts[keep, ]))
  y <- estimateDisp(y, design)
  f <- glmQLFit(y, design)
  p <- rep(NA, length(keep)); p[keep] <- p.adjust(glmQLFTest(f, coef = 2)$table$PValue, "BH")
  out$edgeR <- p
  # limma-voom
  v <- voom(y, design)
  fv <- eBayes(lmFit(v, design))
  p <- rep(NA, length(keep)); p[keep] <- topTable(fv, coef = 2, number = Inf, sort.by = "none")$adj.P.Val
  out$voom <- p
  out
}

score <- function(padj, truth, alpha = .05) {
  call <- !is.na(padj) & padj < alpha
  c(discoveries = sum(call),
    observed_FDR = if (sum(call)) sum(call & !truth) / sum(call) else 0,
    sensitivity = sum(call & truth) / sum(truth))
}

# -----------------------------------------------------------------------------
banner("3. Power and FDR vs replicates (this is the slow part)")
# -----------------------------------------------------------------------------
grid <- expand.grid(n = c(2, 3, 5, 8), rep = 1:2,
                    variability = c("cell line (BCV~0.1)", "tissue (BCV~0.3)"),
                    stringsAsFactors = FALSE)
perf <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  sim <- simulate(grid$n[i], disp_mult = ifelse(grepl("tissue", grid$variability[i]), 9, 1))
  r <- run_methods(sim)
  do.call(rbind, lapply(names(r), function(m)
    data.frame(grid[i, ], method = m, t(score(r[[m]], sim$truth)))))
}))
agg <- aggregate(cbind(observed_FDR, sensitivity) ~ n + method + variability, perf, mean)
print(agg[order(agg$variability, agg$method, agg$n), ], digits = 3, row.names = FALSE)
write.csv(perf, file.path(od, "simulation_performance.csv"), row.names = FALSE)

p <- ggplot(agg, aes(n, sensitivity, colour = method)) + geom_line() + geom_point() +
  facet_wrap(~ variability) +
  labs(x = "replicates per group", y = "sensitivity (true DE found)",
       title = "Power grows with replicates, not with sequencing depth alone")
p2 <- ggplot(agg, aes(n, observed_FDR, colour = method)) + geom_line() + geom_point() +
  facet_wrap(~ variability) +
  geom_hline(yintercept = .05, linetype = 2) +
  labs(x = "replicates per group", y = "observed FDR", title = "Is the nominal 5% FDR honest?")
save_plot(p, file.path(od, "power_curve.png"), width = 9, height = 4)
save_plot(p2, file.path(od, "fdr_curve.png"), width = 9, height = 4)
cat("\nTakeaways: with clean cell lines even n=2-3 finds most changes. With\n",
    "tissue-level variability the same n misses a large share, and only more\n",
    "replicates recover it. Estimate YOUR BCV (module 05) before choosing n.\n",
    "Check the FDR column too: edgeR-QL and voom stay near 5%, while DESeq2's\n",
    "Wald test can run above the nominal FDR with few, noisy replicates (QL methods\n",
    "are deliberately more conservative). Exact numbers move between runs -- change\n",
    "the seed or add reps to see how much. Real data has outliers and batch effects\n",
    "the simulator lacks, so this is the BEST case.\n")

# Sensitivity depends on effect size and expression
sim <- simulate(3, disp_mult = 9); r <- run_methods(sim)
dd <- data.frame(abs_lfc = abs(sim$true_lfc), found = !is.na(r$DESeq2) & r$DESeq2 < .05,
                 truth = sim$truth)[sim$truth, ]
dd$lfc_bin <- cut(dd$abs_lfc, c(0.5, 1, 1.5, 2, 3, Inf), include.lowest = TRUE)
cat("\nDESeq2 sensitivity by true |log2FC| (n=3, tissue-level variability):\n")
print(round(tapply(dd$found, dd$lfc_bin, mean), 2))

# -----------------------------------------------------------------------------
banner("4. Pitfall: technical replicates are not biological replicates")
# -----------------------------------------------------------------------------
# Mock comparison again (2 donors vs 2 donors), but now pretend each library
# was sequenced 3 times and treat each run as a separate sample.
tech <- do.call(cbind, lapply(1:3, function(i)
  matrix(rpois(length(cts_c), cts_c), nrow = nrow(cts_c))))   # Poisson = sequencing noise
colnames(tech) <- paste0(rep(colnames(cts_c), 3), "_run", rep(1:3, each = 4))
rownames(tech) <- rownames(cts_c)
cd_t <- data.frame(fake = rep(cd_c$fake, 3), donor = rep(cd_c$cell, 3), row.names = colnames(tech))
dds_t <- DESeq(DESeqDataSetFromMatrix(tech, cd_t, ~ fake), quiet = TRUE)
n_wrong <- sum(results(dds_t)$padj < .05, na.rm = TRUE)
dds_ok <- collapseReplicates(DESeqDataSetFromMatrix(tech, cd_t, ~ fake), groupby = cd_t$donor)
n_right <- sum(results(DESeq(dds_ok, quiet = TRUE))$padj < .05, na.rm = TRUE)
cat("Mock comparison, technical runs as samples: ", n_wrong, "'DE' genes\n",
    "Same data, runs summed per donor (collapseReplicates):", n_right, "\n",
    "Technical replicates only measure sequencing noise, so the model thinks the\n",
    "data are far more precise than they are. Sum them; never count them as n.\n")

# -----------------------------------------------------------------------------
banner("5. Pitfall: double dipping")
# -----------------------------------------------------------------------------
# Pure noise: 6 samples from the SAME distribution, random labels
sim0 <- simulate(3, frac_de = 0)
d0 <- DESeq(DESeqDataSetFromMatrix(sim0$counts, data.frame(group = sim0$group), ~ group),
            quiet = TRUE)
r0 <- results(d0)
top <- head(order(r0$pvalue), 50)
cat("True DE genes: 0. Smallest padj:", signif(min(r0$padj, na.rm = TRUE), 2), "\n")
m <- assay(vst(d0, blind = TRUE))[top, ]
m <- m - rowMeans(m)
save_plot(function() pheatmap(m, annotation_col = data.frame(group = sim0$group,
                                                             row.names = colnames(m)),
                              show_rownames = FALSE,
                              main = "Top-50 p-value genes from PURE NOISE"),
          file.path(od, "double_dipping_heatmap.png"), width = 6, height = 6)
cat("The heatmap clusters the samples perfectly by group -- from noise. Genes\n",
    "SELECTED for separating the groups will always separate the groups. Such a\n",
    "heatmap illustrates a result; it is never evidence for one.\n")
