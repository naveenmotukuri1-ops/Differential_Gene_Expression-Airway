# =============================================================================
# MODULE 02  --  Quality control and exploratory analysis   (airway, human)
# =============================================================================
# Learning goals
#   * Know what a count matrix is and what to check BEFORE any statistics
#   * Library size, detection, gene biotype composition
#   * Why we filter low counts, and how to choose a filter
#   * Why raw counts can't be used for PCA/clustering: the mean-variance link
#   * log2(x+1) vs VST vs rlog
#   * PCA, sample-to-sample distances, what "the samples cluster by donor" means
#
# Run from the course root:  Rscript 02_qc_exploration.R
# Outputs: results/02_qc/
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
})
od <- out_dir("02_qc")

aw <- load_airway()
cts <- aw$counts; cd <- aw$coldata; genes <- aw$genes

# -----------------------------------------------------------------------------
banner("1. What are we looking at?")
# -----------------------------------------------------------------------------
cat("Genes x samples:", dim(cts), "\n")
print(cd)
cat("\nThe design: every donor (cell) contributes one control AND one treated\n",
    "sample. This is a PAIRED design. Remember it: it decides the model in 04.\n")
print(table(cd$cell, cd$dex))

# Rule #1 of the field: confirm sample order matches between counts and metadata
stopifnot(all(colnames(cts) == rownames(cd)))

# -----------------------------------------------------------------------------
banner("2. Library size and detection")
# -----------------------------------------------------------------------------
qc <- data.frame(sample = colnames(cts),
                 dex = cd$dex, cell = cd$cell,
                 lib_size_M = colSums(cts) / 1e6,
                 genes_detected = colSums(cts > 0),
                 pct_top100 = apply(cts, 2, function(x) 100 * sum(sort(x, TRUE)[1:100]) / sum(x)))
print(qc, digits = 3)
write.csv(qc, file.path(od, "sample_qc.csv"), row.names = FALSE)

p <- ggplot(qc, aes(sample, lib_size_M, fill = dex)) + geom_col() +
  geom_hline(yintercept = mean(qc$lib_size_M), linetype = 2) +
  labs(y = "Assigned reads (millions)", x = NULL, title = "Library sizes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, file.path(od, "library_sizes.png"))
cat("\nLibraries range", round(min(qc$lib_size_M), 1), "-", round(max(qc$lib_size_M), 1),
    "M reads. A 2-fold range is normal; normalization (module 03) handles it.\n",
    "A sample with 10x fewer reads, or far fewer genes detected, is a red flag.\n")

# Biotype composition: rRNA or mitochondrial dominance points to library problems
bt <- rowsum(cts, genes$biotype)
bt_pct <- sweep(bt, 2, colSums(bt), "/") * 100
top_bt <- names(sort(rowMeans(bt_pct), decreasing = TRUE))[1:6]
cat("\n% of counts by gene biotype (top 6):\n")
print(round(bt_pct[top_bt, ], 2))
mt <- colSums(cts[genes$chr == "MT", ]) / colSums(cts) * 100
cat("\n% counts on mitochondrial genes:\n"); print(round(mt, 2))

# -----------------------------------------------------------------------------
banner("3. Filtering low-count genes")
# -----------------------------------------------------------------------------
cat("Genes with zero counts in ALL samples:", sum(rowSums(cts) == 0), "of", nrow(cts), "\n")
# A common, defensible rule: keep genes with >= 10 counts in at least as many
# samples as the smallest group (here 4). edgeR::filterByExpr does a smarter
# version of this based on CPM; we'll meet it in module 05.
keep <- rowSums(cts >= 10) >= 4
cat("Kept by 'count >= 10 in >= 4 samples':", sum(keep), "\n")
cat("Why filter? Genes with almost no reads carry no information about DE,\n",
    "but they DO count towards multiple-testing correction and slow everything down.\n",
    "Note: DESeq2's independent filtering does this automatically at results()\n",
    "time, so this pre-filter is mainly for speed and plots.\n")

log_cpm <- log2(sweep(cts, 2, colSums(cts), "/") * 1e6 + 1)
dens <- data.frame(value = c(log_cpm), sample = rep(colnames(cts), each = nrow(cts)),
                   kept = rep(keep, ncol(cts)))
p <- ggplot(dens, aes(value, colour = sample)) + geom_density() +
  facet_wrap(~ ifelse(kept, "after filter", "removed genes"), scales = "free_y") +
  labs(x = "log2(CPM + 1)", title = "Distribution of expression, filtered vs removed")
save_plot(p, file.path(od, "filtering_density.png"), width = 10)

cts_f <- cts[keep, ]

# -----------------------------------------------------------------------------
banner("4. The mean-variance relationship (why transformations matter)")
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(cts_f, cd, design = ~ cell + dex)
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized = TRUE)

mv <- data.frame(mean = rowMeans(norm), var = rowVars(norm))
p <- ggplot(mv, aes(mean, var)) + geom_point(alpha = .15, size = .5) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Variance grows faster than the mean (red line = Poisson)",
       x = "mean normalized count", y = "variance")
save_plot(p, file.path(od, "mean_variance.png"))
cat("Counts are over-dispersed: variance > mean (the red Poisson line).\n",
    "This is why we model counts with the Negative Binomial distribution.\n")

# Three transformations, and what each does to the low-count end
log2c <- log2(norm + 1)
vsd   <- vst(dds, blind = TRUE)        # blind = TRUE for unsupervised QC
rld   <- rlog(dds, blind = TRUE)       # slow for >30 samples; fine here

sd_plot <- function(m, label) {
  data.frame(rank = rank(rowMeans(m)), sd = rowSds(m), method = label)
}
sdd <- rbind(sd_plot(log2c, "log2(norm+1)"), sd_plot(assay(vsd), "VST"),
             sd_plot(assay(rld), "rlog"))
p <- ggplot(sdd, aes(rank, sd)) + geom_hex(bins = 60) + geom_smooth(colour = "red", se = FALSE) +
  facet_wrap(~ method) + labs(x = "rank of mean expression", y = "SD across samples",
                              title = "Which transform makes SD independent of the mean?")
if (!requireNamespace("hexbin", quietly = TRUE)) p$layers[[1]] <- geom_point(alpha = .1, size = .3)
save_plot(p, file.path(od, "meanSd_transforms.png"), width = 11, height = 4)
cat("log2(x+1) inflates SD of low-count genes (they dominate PCA with noise).\n",
    "VST and rlog flatten that. Use VST by default; rlog when size factors vary a lot.\n")

# -----------------------------------------------------------------------------
banner("5. PCA")
# -----------------------------------------------------------------------------
pca_df <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE, ntop = 500)
pv <- round(100 * attr(pca_df, "percentVar"))
p <- ggplot(pca_df, aes(PC1, PC2, colour = dex, shape = cell)) + geom_point(size = 4) +
  labs(x = paste0("PC1: ", pv[1], "%"), y = paste0("PC2: ", pv[2], "%"),
       title = "PCA of VST data, top 500 variable genes")
save_plot(p, file.path(od, "pca_vst.png"))

# Do it by hand once so plotPCA is not a black box
v <- assay(vsd)
top <- head(order(rowVars(v), decreasing = TRUE), 500)
pc <- prcomp(t(v[top, ]))
cat("Variance explained by PC1-4 (manual prcomp):",
    round(100 * pc$sdev[1:4]^2 / sum(pc$sdev^2), 1), "%\n")
cat("Look at the plot: PC1 separates treatment; PC2 is donor (N080611 sits apart).\n",
    "Within each donor, the treated sample shifts the same way along PC1.\n",
    "Donor is a strong source of variation -> put `cell` in the design (module 04).\n")

# Which genes drive PC1?
load1 <- sort(pc$rotation[, 1])
drivers <- c(head(load1, 5), tail(load1, 5))
cat("\nGenes with the largest PC1 loadings:\n")
print(data.frame(symbol = genes[names(drivers), "symbol"], loading = round(drivers, 3)))

# -----------------------------------------------------------------------------
banner("6. Sample-to-sample distances")
# -----------------------------------------------------------------------------
sd_mat <- as.matrix(dist(t(assay(vsd))))
ann <- data.frame(dex = cd$dex, cell = cd$cell, row.names = colnames(cts))
save_plot(function() pheatmap(sd_mat, annotation_col = ann, annotation_row = ann,
                              col = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
                              main = "Euclidean distance on VST", silent = FALSE),
          file.path(od, "sample_distances.png"), width = 7, height = 6)

# Correlation is a useful (but less sensitive) alternative; values near 1 are expected
cor_mat <- cor(assay(vsd), method = "pearson")
cat("Min pairwise Pearson r between samples:", round(min(cor_mat), 4), "\n")

# -----------------------------------------------------------------------------
banner("7. Save for later modules")
# -----------------------------------------------------------------------------
saveRDS(list(dds = dds, vsd = vsd), file.path(od, "airway_dds_vsd.rds"))
cat("Done. Look at every PNG in", od, "and answer the exercises in GUIDE.md.\n")
