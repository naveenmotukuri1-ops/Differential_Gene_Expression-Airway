# =============================================================================
# MODULE 03  --  Normalization, from first principles   (airway, human)
# =============================================================================
# Learning goals
#   * Why "divide by library size" (CPM) is not enough: composition bias
#   * DESeq2 median-of-ratios: implement it by hand, then check against DESeq2
#   * edgeR TMM: same idea, different estimator
#   * TPM / RPKM: what gene length does, and why they are NOT inputs to DE tools
#   * A simulated example where CPM gives false DE and median-of-ratios doesn't
#
# Run:  Rscript 03_normalization.R        Outputs: results/03_norm/
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
})
od <- out_dir("03_norm")

aw <- load_airway()
keep <- rowSums(aw$counts >= 10) >= 4
cts <- aw$counts[keep, ]; cd <- aw$coldata; genes <- aw$genes[keep, ]

# -----------------------------------------------------------------------------
banner("1. Median-of-ratios by hand (DESeq2's method)")
# -----------------------------------------------------------------------------
# Step 1: geometric mean of each gene across samples = a pseudo-reference sample
log_geo_means <- rowMeans(log(cts))           # -Inf if any zero -> gene dropped
use <- is.finite(log_geo_means)
cat("Genes usable (no zeros in any sample):", sum(use), "\n")
# Step 2: ratio of each sample to the pseudo-reference, per gene
log_ratios <- log(cts[use, ]) - log_geo_means[use]
# Step 3: size factor = median ratio. The median is robust to the minority of
# genes that are truly DE -- that is the whole trick.
sf_manual <- exp(apply(log_ratios, 2, median))

dds <- DESeqDataSetFromMatrix(cts, cd, ~ cell + dex)
dds <- estimateSizeFactors(dds)
cmp <- data.frame(manual = sf_manual, DESeq2 = sizeFactors(dds),
                  lib_size_ratio = colSums(cts) / exp(mean(log(colSums(cts)))))
print(round(cmp, 4))
stopifnot(isTRUE(all.equal(unname(sf_manual), unname(sizeFactors(dds)))))
cat("Manual size factors match DESeq2 exactly.\n")

# Visualize the ratio distribution for two samples
rat <- data.frame(log2ratio = c(log_ratios[, 1], log_ratios[, 6]) / log(2),
                  sample = rep(colnames(cts)[c(1, 6)], each = sum(use)))
p <- ggplot(rat, aes(log2ratio)) + geom_histogram(bins = 100) +
  facet_wrap(~ sample) + geom_vline(data = data.frame(sample = colnames(cts)[c(1, 6)],
                                                      m = log2(sf_manual[c(1, 6)])),
                                    aes(xintercept = m), colour = "red") +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(title = "Per-gene log2 ratio to pseudo-reference; red = size factor")
save_plot(p, file.path(od, "median_of_ratios.png"), width = 9, height = 4)

# -----------------------------------------------------------------------------
banner("2. TMM (edgeR's method)")
# -----------------------------------------------------------------------------
y <- DGEList(cts, group = cd$dex)
y <- calcNormFactors(y, method = "TMM")
# edgeR's "effective library size" = lib.size * norm.factor; scale to compare
eff <- y$samples$lib.size * y$samples$norm.factors
cmp$TMM_effective <- eff / exp(mean(log(eff)))
print(round(cmp, 4))
cat("Correlation DESeq2 vs TMM scaling:", round(cor(cmp$DESeq2, cmp$TMM_effective), 4), "\n",
    "The two methods agree closely on well-behaved data. They differ when a few\n",
    "genes take a big share of reads, or when composition shifts globally.\n")

# -----------------------------------------------------------------------------
banner("3. Composition bias: a simulation you should never forget")
# -----------------------------------------------------------------------------
# Take sample 1, make a copy, and in the copy massively up-regulate 5% of genes
# (like a highly expressed transgene or a contamination). Nothing else changes.
set.seed(1)
base <- cts[, 1]
sim <- cbind(A1 = rpois(length(base), base), A2 = rpois(length(base), base),
             B1 = rpois(length(base), base), B2 = rpois(length(base), base))
up <- sample(nrow(sim), round(0.05 * nrow(sim)))
sim[up, c("B1", "B2")] <- sim[up, c("B1", "B2")] * 10
truth_null <- setdiff(seq_len(nrow(sim)), up)     # 95% genes are truly unchanged

cpm_sim <- cpm(sim)                              # library size only
dds_sim <- estimateSizeFactors(DESeqDataSetFromMatrix(sim,
                                 data.frame(g = factor(c("A", "A", "B", "B"))), ~ g))
mor_sim <- counts(dds_sim, normalized = TRUE)

lfc_cpm <- log2(rowMeans(cpm_sim[, 3:4]) + 1) - log2(rowMeans(cpm_sim[, 1:2]) + 1)
lfc_mor <- log2(rowMeans(mor_sim[, 3:4]) + 1) - log2(rowMeans(mor_sim[, 1:2]) + 1)
cat("Median log2FC of the 95% UNCHANGED genes:\n",
    "  CPM (total-count) normalization :", round(median(lfc_cpm[truth_null]), 3), "\n",
    "  median-of-ratios               :", round(median(lfc_mor[truth_null]), 3), "\n",
    "CPM makes every unchanged gene look DOWN-regulated, because the up-regulated\n",
    "genes 'stole' sequencing depth. Median-of-ratios (and TMM) assume most genes\n",
    "are not DE and correct for it.\n")
df <- data.frame(lfc = c(lfc_cpm[truth_null], lfc_mor[truth_null]),
                 method = rep(c("CPM", "median-of-ratios"), each = length(truth_null)))
p <- ggplot(df, aes(lfc, fill = method)) + geom_density(alpha = .5) +
  geom_vline(xintercept = 0) + coord_cartesian(xlim = c(-2, 1.5)) +
  labs(x = "log2 fold change of truly unchanged genes", title = "Composition bias")
save_plot(p, file.path(od, "composition_bias_sim.png"))

# -----------------------------------------------------------------------------
banner("4. Gene length: RPKM/FPKM and TPM")
# -----------------------------------------------------------------------------
len_kb <- genes$exonic_length / 1000
rpk  <- cts / len_kb
tpm  <- sweep(rpk, 2, colSums(rpk), "/") * 1e6
rpkm <- sweep(cts, 2, colSums(cts) / 1e6, "/") / len_kb
cat("Column sums -- TPM always sums to 1e6, RPKM does not:\n")
print(round(rbind(TPM = colSums(tpm), RPKM = colSums(rpkm)), 0))
cat("\nTPM is good for comparing genes WITHIN a sample (gene A vs gene B).\n",
    "For DE BETWEEN samples, gene length cancels out (same gene, same length), so\n",
    "DE tools take RAW counts + size factors. Never give TPM/FPKM to DESeq2/edgeR:\n",
    "they need the counts to model the count noise.\n")

# Longer genes collect more reads at the same expression level:
lc <- data.frame(len = genes$exonic_length, mean_count = rowMeans(cts),
                 mean_tpm = rowMeans(tpm))
cat("Spearman cor(length, raw count) =", round(cor(lc$len, lc$mean_count, method = "s"), 3),
    "  cor(length, TPM) =", round(cor(lc$len, lc$mean_tpm, method = "s"), 3), "\n")

# -----------------------------------------------------------------------------
banner("5. Which normalized values for what?")
# -----------------------------------------------------------------------------
cat("  DE testing ...................... raw counts (tool normalizes internally)\n",
    " PCA / clustering / heatmaps ..... VST or rlog (DESeq2), logCPM w/ prior (edgeR)\n",
    " Plotting one gene ............... normalized counts, log scale\n",
    " Comparing genes in one sample ... TPM\n",
    " Sharing with collaborators ...... raw counts + metadata (+ TPM for browsing)\n")

norm_tab <- data.frame(gene_id = rownames(cts), symbol = genes$symbol,
                       round(counts(dds, normalized = TRUE), 1), check.names = FALSE)
write.csv(norm_tab, file.path(od, "normalized_counts_DESeq2.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = rownames(cts), symbol = genes$symbol, round(tpm, 2),
                     check.names = FALSE), file.path(od, "TPM.csv"), row.names = FALSE)
cat("Wrote normalized counts and TPM tables to", od, "\n")
