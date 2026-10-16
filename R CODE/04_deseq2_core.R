# =============================================================================
# MODULE 04  --  Differential expression with DESeq2   (airway, human)
# =============================================================================
# Learning goals
#   * Design formulas: why ~ cell + dex and not ~ dex (paired/blocked design)
#   * What DESeq() does: size factors -> dispersions -> GLM -> Wald test
#   * Dispersion shrinkage: read the dispersion plot
#   * results(): contrasts, alpha, independent filtering, Cook's outliers
#   * log2FC shrinkage (apeglm) and why you rank/plot with it
#   * Testing against a fold-change threshold (lfcThreshold) vs post-hoc cutoffs
#   * MA plot, volcano, heatmap of top genes, single-gene plots
#
# Run:  Rscript 04_deseq2_core.R          Outputs: results/04_deseq2/
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(ggrepel)
})
od <- out_dir("04_deseq2")

aw <- load_airway()
keep <- rowSums(aw$counts >= 10) >= 4
cts <- aw$counts[keep, ]; cd <- aw$coldata; genes <- aw$genes[keep, ]

# -----------------------------------------------------------------------------
banner("1. Design: the single most important decision")
# -----------------------------------------------------------------------------
# The LAST term is what results() tests by default. Blocking factors go first.
dds_naive  <- DESeqDataSetFromMatrix(cts, cd, design = ~ dex)
dds        <- DESeqDataSetFromMatrix(cts, cd, design = ~ cell + dex)
cat("Reference level of dex:", levels(dds$dex)[1],
    " (log2FC will be treated / control)\n")

# -----------------------------------------------------------------------------
banner("2. Run DESeq()")
# -----------------------------------------------------------------------------
dds       <- DESeq(dds)
dds_naive <- DESeq(dds_naive, quiet = TRUE)
print(resultsNames(dds))

save_plot(function() plotDispEsts(dds, main = "Dispersion estimates"),
          file.path(od, "dispersion.png"), width = 6, height = 5)
cat("\nDispersion plot: black = per-gene estimate, red = fitted trend,\n",
    "blue = shrunken (final) values. Circled = outliers kept at their own value.\n",
    "Shrinkage borrows information across genes -- crucial with n=4 per group.\n")

# -----------------------------------------------------------------------------
banner("3. Results: paired vs unpaired design")
# -----------------------------------------------------------------------------
res       <- results(dds, contrast = c("dex", "treated", "control"), alpha = 0.05)
res_naive <- results(dds_naive, contrast = c("dex", "treated", "control"), alpha = 0.05)
summary(res)
cat("Significant at FDR 5%:  paired design =", sum(res$padj < .05, na.rm = TRUE),
    "  unpaired design =", sum(res_naive$padj < .05, na.rm = TRUE), "\n",
    "Modelling donor removes donor-to-donor variation from the residual noise,\n",
    "so the same data gives more power. Ignoring structure = wasting samples.\n")

# What the columns mean
cat("\nColumns of the results table:\n"); print(mcols(res)$description)

# -----------------------------------------------------------------------------
banner("4. Independent filtering and outliers")
# -----------------------------------------------------------------------------
md <- metadata(res)
cat("Independent filtering threshold (baseMean):", round(md$filterThreshold, 2), "\n",
    "Genes with NA padj because of low counts:",
    sum(!is.na(res$pvalue) & is.na(res$padj)), "\n",
    "Genes with NA p-value (Cook's outlier or all zero):", sum(is.na(res$pvalue)), "\n")
save_plot(function() {
  plot(md$filterNumRej, type = "b", xlab = "quantile of baseMean filtered",
       ylab = "number of rejections", main = "Independent filtering")
  lines(md$lo.fit, col = "red"); abline(v = md$filterTheta)
}, file.path(od, "independent_filtering.png"), width = 6, height = 5)

# p-value histogram: the first diagnostic I look at for ANY DE analysis
save_plot(ggplot(as.data.frame(res[res$baseMean > 1, ]), aes(pvalue)) +
            geom_histogram(breaks = seq(0, 1, .025)) +
            labs(title = "p-value histogram",
                 subtitle = "Flat + spike near 0 = healthy. U-shape / hump = model problem."),
          file.path(od, "pvalue_hist.png"))

# -----------------------------------------------------------------------------
banner("5. log2 fold-change shrinkage")
# -----------------------------------------------------------------------------
# Low-count genes have noisy, inflated LFCs. apeglm shrinks them towards 0
# in proportion to their uncertainty; well-measured genes barely move.
coef_name <- "dex_treated_vs_control"
res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)

save_plot(function() {
  par(mfrow = c(1, 2))
  plotMA(res, ylim = c(-5, 5), main = "MLE log2FC")
  plotMA(res_shr, ylim = c(-5, 5), main = "apeglm-shrunken log2FC")
}, file.path(od, "MA_mle_vs_shrunken.png"), width = 11, height = 5)

lf <- data.frame(mle = res$log2FoldChange, shr = res_shr$log2FoldChange,
                 baseMean = res$baseMean)
p <- ggplot(lf, aes(mle, shr, colour = log10(baseMean + 1))) + geom_point(size = .5) +
  geom_abline(linetype = 2) + scale_colour_viridis_c() +
  labs(x = "MLE log2FC", y = "apeglm log2FC", title = "Shrinkage acts mostly on low-count genes")
save_plot(p, file.path(od, "lfc_shrinkage_scatter.png"))
cat("Use shrunken LFCs for ranking, plots, and GSEA. Use the Wald p-values from\n",
    "results() for significance (apeglm s-values are an alternative).\n")

# -----------------------------------------------------------------------------
banner("6. Testing a biologically meaningful effect size")
# -----------------------------------------------------------------------------
# Filtering padj<0.05 & |LFC|>1 AFTER the test does NOT control FDR for the
# claim "|LFC|>1". Put the threshold INTO the null hypothesis instead:
res_lfc1 <- results(dds, contrast = c("dex", "treated", "control"),
                    lfcThreshold = 1, altHypothesis = "greaterAbs", alpha = 0.05)
n_post <- sum(res$padj < .05 & abs(res$log2FoldChange) > 1, na.rm = TRUE)
n_test <- sum(res_lfc1$padj < .05, na.rm = TRUE)
cat("padj<.05 then |LFC|>1 post-hoc:", n_post, "genes\n",
    "Formal test H0: |LFC| <= 1     :", n_test, "genes\n",
    "The formal test is stricter -- and it is the one whose FDR is valid.\n")

# -----------------------------------------------------------------------------
banner("7. The results table you would hand to a biologist")
# -----------------------------------------------------------------------------
tab <- data.frame(gene_id = rownames(res), symbol = genes[rownames(res), "symbol"],
                  biotype = genes[rownames(res), "biotype"],
                  baseMean = round(res$baseMean, 1),
                  log2FC_shrunk = round(res_shr$log2FoldChange, 3),
                  log2FC_mle = round(res$log2FoldChange, 3),
                  lfcSE = round(res$lfcSE, 3), pvalue = res$pvalue, padj = res$padj,
                  padj_lfc1 = res_lfc1$padj)
tab <- tab[order(tab$padj), ]
write.csv(tab, file.path(od, "DE_dex_treated_vs_control.csv"), row.names = FALSE)
cat("Top 15 genes:\n")
print(head(tab[, c("symbol", "baseMean", "log2FC_shrunk", "padj")], 15), row.names = FALSE)
cat("\nSanity check with biology: glucocorticoid targets such as FKBP5, TSC22D3,\n",
    "PER1, ZBTB16, KLF15 and DUSP1 should be up. If they aren't, suspect a sample swap.\n")
print(tab[tab$symbol %in% c("FKBP5", "TSC22D3", "PER1", "ZBTB16", "KLF15", "DUSP1"),
          c("symbol", "log2FC_shrunk", "padj")], row.names = FALSE)

# -----------------------------------------------------------------------------
banner("8. Volcano, heatmap, single genes")
# -----------------------------------------------------------------------------
vt <- tab[!is.na(tab$padj), ]
vt$sig <- ifelse(vt$padj < .05 & abs(vt$log2FC_shrunk) > 1, "FDR<5% & |LFC|>1", "n.s.")
lab <- head(vt[vt$sig != "n.s.", ], 20)
p <- ggplot(vt, aes(log2FC_shrunk, -log10(pvalue), colour = sig)) +
  geom_point(size = .6, alpha = .6) +
  scale_colour_manual(values = c("firebrick", "grey60")) +
  geom_text_repel(data = lab, aes(label = symbol), size = 3, colour = "black", max.overlaps = 30) +
  labs(x = "shrunken log2 fold change", title = "Dexamethasone vs control")
save_plot(p, file.path(od, "volcano.png"), width = 8, height = 6)

vsd <- vst(dds, blind = FALSE)   # blind=FALSE: use the design when we're past QC
top <- head(tab$gene_id, 40)
mat <- assay(vsd)[top, ]
mat <- mat - rowMeans(mat)       # centre each gene: show relative change
rownames(mat) <- genes[top, "symbol"]
ann <- data.frame(dex = cd$dex, cell = cd$cell, row.names = colnames(mat))
save_plot(function() pheatmap(mat, annotation_col = ann, show_rownames = TRUE, fontsize_row = 7,
                              main = "Top 40 DE genes (VST, row-centred)"),
          file.path(od, "heatmap_top40.png"), width = 7, height = 8)

g <- tab$gene_id[1]
gd <- plotCounts(dds, g, intgroup = c("dex", "cell"), returnData = TRUE)
p <- ggplot(gd, aes(dex, count, colour = cell, group = cell)) + geom_point(size = 3) +
  geom_line() + scale_y_log10() +
  labs(title = paste("Top gene:", genes[g, "symbol"]), y = "normalized count",
       subtitle = "Lines join the two samples of each donor = the pairing the model uses")
save_plot(p, file.path(od, "top_gene_paired.png"))

saveRDS(list(dds = dds, res = res, res_shr = res_shr, tab = tab),
        file.path(od, "airway_deseq2.rds"))
cat("Saved DESeq2 objects for modules 05-06 and 08.\n")
