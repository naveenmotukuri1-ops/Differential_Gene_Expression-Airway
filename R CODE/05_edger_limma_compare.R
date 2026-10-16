# =============================================================================
# MODULE 05  --  edgeR quasi-likelihood and limma-voom, compared to DESeq2
# =============================================================================
# Learning goals
#   * The edgeR workflow: DGEList, filterByExpr, TMM, model.matrix, QL F-test
#   * The limma-voom workflow: precision weights turn counts into "microarray-like"
#     data; why voom is fast and handles complex designs/random effects
#   * Design matrices by hand (model.matrix) and makeContrasts
#   * duplicateCorrelation: donor as a RANDOM effect instead of a fixed block
#   * How much do the three methods agree, and where do they disagree?
#
# Run:  Rscript 05_edger_limma_compare.R  Outputs: results/05_compare/
# Needs results from 04 (run 04_deseq2_core.R first).
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
})
od <- out_dir("05_compare")

aw <- load_airway()
cts <- aw$counts; cd <- aw$coldata; genes <- aw$genes
f04 <- file.path(RESULTS_DIR, "04_deseq2", "airway_deseq2.rds")
if (!file.exists(f04)) stop("Run 04_deseq2_core.R first")
de2 <- readRDS(f04)

# -----------------------------------------------------------------------------
banner("1. Design matrix by hand")
# -----------------------------------------------------------------------------
design <- model.matrix(~ cell + dex, data = cd)
print(design)
cat("Columns: intercept = control in reference donor; cellN0xxxx = donor offsets;\n",
    "dextreated = the treatment effect WITHIN donor -- the coefficient we test.\n")

# -----------------------------------------------------------------------------
banner("2. edgeR quasi-likelihood pipeline")
# -----------------------------------------------------------------------------
y <- DGEList(cts, samples = cd, genes = genes)
keep <- filterByExpr(y, design)          # CPM-based, adapts to library size and group size
cat("filterByExpr keeps", sum(keep), "genes\n")
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)                     # TMM (calcNormFactors in older edgeR)
y <- estimateDisp(y, design)
cat("Common BCV (biological coefficient of variation):", round(sqrt(y$common.dispersion), 3),
    "\n  ~0.1-0.2 for cell lines, 0.2-0.4 for inbred animals, 0.4+ for humans\n")
save_plot(function() plotBCV(y), file.path(od, "edgeR_BCV.png"), width = 6, height = 5)

fit <- glmQLFit(y, design, robust = TRUE)
save_plot(function() plotQLDisp(fit), file.path(od, "edgeR_QLdisp.png"), width = 6, height = 5)
qlf <- glmQLFTest(fit, coef = "dextreated")
tt_edger <- topTags(qlf, n = Inf)$table
cat("edgeR QL: FDR<5% =", sum(tt_edger$FDR < .05), "\n")
summary(decideTests(qlf))

# glmTreat = edgeR's version of DESeq2's lfcThreshold
tr <- glmTreat(fit, coef = "dextreated", lfc = 1)
cat("glmTreat |LFC|>1 at FDR 5%:", sum(p.adjust(tr$table$PValue, "BH") < .05), "\n")

# -----------------------------------------------------------------------------
banner("3. limma-voom")
# -----------------------------------------------------------------------------
v <- voom(y, design, plot = FALSE, save.plot = TRUE)
save_plot(function() {
  plot(v$voom.xy, pch = 16, cex = .25, xlab = v$voom.xy$xlab, ylab = v$voom.xy$ylab,
       main = "voom: mean-variance trend")
  lines(v$voom.line, col = "red")
}, file.path(od, "voom_trend.png"), width = 6, height = 5)
cat("voom estimates the mean-variance trend on log-CPM and turns it into a\n",
    "precision weight per observation. Then ordinary limma (with empirical Bayes)\n",
    "applies -- the whole limma toolbox (contrasts, random effects, camera) opens up.\n")
fitv <- lmFit(v, design)
fitv <- eBayes(fitv, robust = TRUE)
tt_voom <- topTable(fitv, coef = "dextreated", number = Inf, sort.by = "none")
cat("limma-voom: FDR<5% =", sum(tt_voom$adj.P.Val < .05), "\n")
save_plot(function() plotSA(fitv, main = "Residual SD vs average log expression"),
          file.path(od, "voom_SA.png"), width = 6, height = 5)

# -----------------------------------------------------------------------------
banner("4. Donor as a random effect (duplicateCorrelation)")
# -----------------------------------------------------------------------------
# Instead of spending 3 df on donor dummies, estimate the within-donor correlation.
# With 4 donors it makes little difference; with many donors/patients and
# unbalanced designs it can matter a lot (and it's needed when donor is
# NOT nested within the comparison, e.g. some donors only in one group).
design_rnd <- model.matrix(~ dex, data = cd)
v2   <- voom(y, design_rnd)
corf <- duplicateCorrelation(v2, design_rnd, block = cd$cell)
cat("Consensus within-donor correlation:", round(corf$consensus, 3), "\n")
v2   <- voom(y, design_rnd, block = cd$cell, correlation = corf$consensus)
corf <- duplicateCorrelation(v2, design_rnd, block = cd$cell)  # iterate once (recommended)
fit_r <- eBayes(lmFit(v2, design_rnd, block = cd$cell, correlation = corf$consensus))
tt_rnd <- topTable(fit_r, coef = "dextreated", number = Inf, sort.by = "none")
cat("voom + random donor effect: FDR<5% =", sum(tt_rnd$adj.P.Val < .05), "\n")

# -----------------------------------------------------------------------------
banner("5. Agreement between DESeq2, edgeR and limma-voom")
# -----------------------------------------------------------------------------
common <- Reduce(intersect, list(rownames(tt_edger), rownames(tt_voom),
                                 rownames(de2$res)[!is.na(de2$res$padj)]))
cmp <- data.frame(gene_id = common, symbol = genes[common, "symbol"],
                  lfc_deseq2 = de2$res[common, "log2FoldChange"],
                  lfc_edger  = tt_edger[common, "logFC"],
                  lfc_voom   = tt_voom[common, "logFC"],
                  padj_deseq2 = de2$res[common, "padj"],
                  fdr_edger   = tt_edger[common, "FDR"],
                  fdr_voom    = tt_voom[common, "adj.P.Val"])
sig <- with(cmp, data.frame(DESeq2 = padj_deseq2 < .05, edgeR = fdr_edger < .05,
                            voom = fdr_voom < .05))
cat("Genes tested by all three:", length(common), "\n")
print(colSums(sig))
cat("Called by all three:", sum(rowSums(sig) == 3),
    "  by exactly one:", sum(rowSums(sig) == 1), "\n")
save_plot(function() vennDiagram(vennCounts(sig), circle.col = c("red", "blue", "darkgreen"),
                                 main = "FDR < 5%"),
          file.path(od, "venn_methods.png"), width = 6, height = 6)
cat("LFC correlation (Pearson):\n")
print(round(cor(cmp[, c("lfc_deseq2", "lfc_edger", "lfc_voom")]), 3))

p <- ggplot(cmp, aes(lfc_deseq2, lfc_edger, colour = rowSums(sig) == 3)) +
  geom_point(size = .5) + geom_abline(linetype = 2) +
  scale_colour_manual(values = c("grey60", "firebrick"), name = "sig in all 3") +
  labs(x = "DESeq2 log2FC (MLE)", y = "edgeR log2FC")
save_plot(p, file.path(od, "lfc_deseq2_vs_edger.png"))

# Where do they disagree? Usually low counts, or one-sample outliers.
only1 <- cmp[rowSums(sig) == 1, ]
only1$baseMean <- de2$res[only1$gene_id, "baseMean"]
cat("\nMedian baseMean -- genes called by only one method:", round(median(only1$baseMean)),
    " vs all-three genes:",
    round(median(de2$res[cmp$gene_id[rowSums(sig) == 3], "baseMean"])), "\n",
    "Disagreements are lower-count genes sitting near the FDR cut-off (check their\n",
    "padj in the methods that missed them: mostly 0.05-0.15). Report the robust core, and\n",
    "never 'method-shop' for the tool that gives the answer you hoped for.\n")
write.csv(cmp, file.path(od, "three_method_comparison.csv"), row.names = FALSE)
saveRDS(list(y = y, v = v, fitv = fitv, design = design), file.path(od, "airway_limma.rds"))
