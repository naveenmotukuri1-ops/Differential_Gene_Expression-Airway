# =============================================================================
# MODULE 07  --  Multi-factor designs, interactions, LRT, hidden batch effects
#                (GSE96870 mouse cerebellum, influenza time course, 22 samples)
# =============================================================================
# Learning goals
#   * Exploring a bigger design: which factor dominates the PCA?
#   * Confounding: why ~ infection + time cannot be fitted here
#   * Multi-level factors and explicit contrasts (Day4 vs Day0, Day8 vs Day4...)
#   * Interactions: "is the infection response different in males and females?"
#   * The "group" trick as an easier-to-read alternative to interaction terms
#   * Likelihood ratio test (LRT): any change over time, like an ANOVA
#   * Surrogate variable analysis (sva) for unknown batch effects
#   * removeBatchEffect: for PLOTS only, never before testing
#
# Run:  Rscript 07_multifactor_designs.R   Outputs: results/07_multifactor/
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(sva)
  library(pheatmap)
})
od <- out_dir("07_multifactor")

ms <- load_mouse()
cts <- ms$counts; cd <- ms$coldata; genes <- ms$genes

# -----------------------------------------------------------------------------
banner("1. The design")
# -----------------------------------------------------------------------------
print(with(cd, table(time, sex)))
print(with(cd, table(time, infection)))
cat("Day0 == NonInfected and Day4/Day8 == InfluenzaA: infection is completely\n",
    "determined by time. Remember this.\n")

keep <- rowSums(cts >= 10) >= 7          # 7 = size of the smallest time group
cts <- cts[keep, ]; genes <- genes[keep, ]
cat("Genes kept:", nrow(cts), "\n")

# -----------------------------------------------------------------------------
banner("2. Exploration: what drives variation?")
# -----------------------------------------------------------------------------
dds0 <- DESeqDataSetFromMatrix(cts, cd, ~ 1)
vsd0 <- vst(dds0, blind = TRUE)
pca <- plotPCA(vsd0, intgroup = c("time", "sex"), returnData = TRUE, ntop = 500)
pv <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = time, shape = sex)) + geom_point(size = 3.5) +
  labs(x = paste0("PC1 ", pv[1], "%"), y = paste0("PC2 ", pv[2], "%"),
       title = "Mouse cerebellum: PCA (VST, top 500)")
save_plot(p, file.path(od, "pca_time_sex.png"))

# Sex check -- a QC step that catches sample swaps in every mammalian dataset
sexg <- intersect(c("Xist", "Ddx3y", "Kdm5d", "Eif2s3y", "Uty"), rownames(vsd0))
sx <- sapply(levels(cd$sex), function(s) rowMeans(assay(vsd0)[sexg, cd$sex == s, drop = FALSE]))
cat("Mean VST of sex-linked genes by annotated sex:\n")
print(round(sx, 2))
cat("Xist high in females, Y genes high in males -> annotation is consistent.\n")

# Variance explained by each factor on top-variable genes (quick ANOVA-style R2)
v <- assay(vsd0); top <- head(order(apply(v, 1, var), decreasing = TRUE), 1000)
r2 <- sapply(c("time", "sex"), function(f) {
  mean(apply(v[top, ], 1, function(g) summary(lm(g ~ cd[[f]]))$r.squared))
})
cat("Mean R^2 over top-1000 variable genes:\n"); print(round(r2, 3))

# -----------------------------------------------------------------------------
banner("3. Confounding: a design that can't be fitted")
# -----------------------------------------------------------------------------
mm <- model.matrix(~ infection + time, cd)
cat("Columns:", ncol(mm), " rank:", qr(mm)$rank, "-> NOT full rank\n")
err <- tryCatch(DESeqDataSetFromMatrix(cts, cd, ~ infection + time),
                error = function(e) conditionMessage(e))
cat("DESeq2 says:\n  ", substr(as.character(err), 1, 200), "...\n")
cat("No statistics can separate two factors that always move together.\n",
    "Solve this at the DESIGN stage of an experiment, not the analysis stage.\n")

# -----------------------------------------------------------------------------
banner("4. Additive model ~ sex + time, with contrasts")
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(cts, cd, ~ sex + time)
dds <- DESeq(dds, quiet = TRUE)
print(resultsNames(dds))
cons <- list(Day4_vs_Day0 = c("time", "Day4", "Day0"),
             Day8_vs_Day0 = c("time", "Day8", "Day0"),
             Day8_vs_Day4 = c("time", "Day8", "Day4"),   # not a coefficient: needs contrast
             Male_vs_Female = c("sex", "Male", "Female"))
res_list <- lapply(cons, function(cn) results(dds, contrast = cn, alpha = .05))
n_sig <- sapply(res_list, function(r) c(up = sum(r$padj < .05 & r$log2FoldChange > 0, na.rm = TRUE),
                                        down = sum(r$padj < .05 & r$log2FoldChange < 0, na.rm = TRUE)))
print(n_sig)
for (nm in names(res_list)) {
  r <- res_list[[nm]]; r <- r[order(r$padj), ]
  write.csv(data.frame(gene = rownames(r), as.data.frame(r)),
            file.path(od, paste0("DE_", nm, ".csv")), row.names = FALSE)
}
r8 <- res_list$Day8_vs_Day0; r8 <- r8[order(r8$padj), ]
cat("\nTop Day8 vs Day0 genes:\n")
print(head(data.frame(gene = rownames(r8), lfc = round(r8$log2FoldChange, 2),
                      padj = signif(r8$padj, 3)), 15), row.names = FALSE)

# Biology check: two programs with different timing
isg    <- c("Irf7", "Isg15", "Ifit1", "Oasl2")      # interferon-stimulated genes
myelin <- c("Plp1", "Mog", "Mbp")                    # oligodendrocyte / myelin genes
bio <- sapply(c("Day4_vs_Day0", "Day8_vs_Day0"), function(nm)
  round(res_list[[nm]][c(isg, myelin), "log2FoldChange"], 2))
rownames(bio) <- c(isg, myelin)
cat("\nlog2FC of interferon-stimulated vs myelin genes:\n"); print(bio)
cat("The interferon response is strongest at Day4 (Irf7/Isg15/Ifit1 near the top\n",
    "of the Day4 list) and fades by Day8, when myelin genes are clearly down.\n",
    "A time course lets you see the ORDER of events -- one snapshot would not.\n")

# Numeric contrast vectors: the general way (works for any linear combination)
rn <- resultsNames(dds)
cvec <- setNames(numeric(length(rn)), rn)
cvec["time_Day8_vs_Day0"] <- 1; cvec["time_Day4_vs_Day0"] <- -1
r_num <- results(dds, contrast = cvec)
stopifnot(isTRUE(all.equal(r_num$log2FoldChange, res_list$Day8_vs_Day4$log2FoldChange)))
cat("\nNumeric contrast Day8 - Day4 gives identical LFCs to the named contrast.\n")

# -----------------------------------------------------------------------------
banner("5. Interaction: does the infection response differ by sex?")
# -----------------------------------------------------------------------------
dds_int <- DESeqDataSetFromMatrix(cts, cd, ~ sex + time + sex:time)
dds_int <- DESeq(dds_int, test = "LRT", reduced = ~ sex + time, quiet = TRUE)
res_int <- results(dds_int, alpha = .05)
cat("Genes with a significant sex x time interaction (LRT, FDR 5%):",
    sum(res_int$padj < .05, na.rm = TRUE), "\n")
print(resultsNames(dds_int))
cat("Interpreting interaction coefficients:\n",
    " time_Day8_vs_Day0          = Day8 effect in FEMALES (reference sex)\n",
    " sexMale.timeDay8           = how much the Day8 effect DIFFERS in males\n",
    " Day8 effect in males       = time_Day8_vs_Day0 + sexMale.timeDay8\n")
ri <- res_int[order(res_int$pvalue), ]
top_int <- head(rownames(ri), 4)
pd <- do.call(rbind, lapply(top_int, function(g) {
  d <- plotCounts(dds_int, g, intgroup = c("time", "sex"), returnData = TRUE); d$gene <- g; d }))
p <- ggplot(pd, aes(time, count, colour = sex, group = sex)) +
  geom_point(position = position_dodge(.3)) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(.3)) +
  scale_y_log10() + facet_wrap(~ gene, scales = "free_y") +
  labs(title = "Strongest sex x time interaction candidates",
       subtitle = "Non-parallel lines = interaction")
save_plot(p, file.path(od, "interaction_top_genes.png"), width = 8, height = 6)
cat("With 3-4 mice per sex x time cell, interaction tests have LOW power.\n",
    "A non-significant interaction is NOT evidence that males and females respond\n",
    "the same. Plan ~2-4x more replicates if interactions are the main question.\n")

# The group trick: same model, easier contrasts
cd$group <- factor(paste(cd$sex, cd$time, sep = "_"))
dds_g <- DESeqDataSetFromMatrix(cts, cd, ~ group)
dds_g <- DESeq(dds_g, quiet = TRUE)
r_m8 <- results(dds_g, contrast = c("group", "Male_Day8", "Male_Day0"))
r_f8 <- results(dds_g, contrast = c("group", "Female_Day8", "Female_Day0"))
cat("\nGroup design: Day8 response in females:", sum(r_f8$padj < .05, na.rm = TRUE),
    "genes, in males:", sum(r_m8$padj < .05, na.rm = TRUE), "genes.\n",
    "Classic trap: 'genes significant in males but not females' is NOT a test of\n",
    "difference. Only the interaction test answers that question.\n")

# -----------------------------------------------------------------------------
banner("6. LRT: any change across the time course")
# -----------------------------------------------------------------------------
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ sex, quiet = TRUE)
res_lrt <- results(dds_lrt, alpha = .05)
cat("Genes changing at ANY time point (LRT, FDR 5%):", sum(res_lrt$padj < .05, na.rm = TRUE), "\n",
    "(The log2FC column of an LRT result is just one coefficient; ignore it --\n",
    " the p-value refers to all time coefficients together.)\n")
saveRDS(list(dds = dds, dds_lrt = dds_lrt, res_lrt = res_lrt, coldata = cd),
        file.path(od, "mouse_dds.rds"))

# -----------------------------------------------------------------------------
banner("7. Unknown batch effects: surrogate variable analysis")
# -----------------------------------------------------------------------------
norm <- counts(dds, normalized = TRUE)
norm <- norm[rowMeans(norm) > 5, ]
mod  <- model.matrix(~ sex + time, cd)
mod0 <- model.matrix(~ sex, cd)
n_sv <- num.sv(log2(norm + 1), mod, method = "be")
cat("Estimated number of surrogate variables:", n_sv, "\n")
n_sv <- max(1, min(n_sv, 2))      # keep it small for a teaching example
svobj <- svaseq(norm, mod, mod0, n.sv = n_sv)
sv <- svobj$sv; colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
cat("Do the SVs track anything we know? Correlation with sample metadata:\n")
print(round(cor(sv, cbind(mouse = cd$mouse, male = as.numeric(cd$sex == "Male"),
                          day = as.numeric(sub("Day", "", cd$time)))), 2))

cd_sv <- cbind(cd, sv)
fml <- as.formula(paste("~", paste(colnames(sv), collapse = " + "), "+ sex + time"))
dds_sv <- DESeqDataSetFromMatrix(cts, cd_sv, fml)
dds_sv <- DESeq(dds_sv, quiet = TRUE)
r_sv <- results(dds_sv, contrast = c("time", "Day8", "Day0"), alpha = .05)
cat("Day8 vs Day0 significant: without SVs =", sum(res_list$Day8_vs_Day0$padj < .05, na.rm = TRUE),
    "  with SVs =", sum(r_sv$padj < .05, na.rm = TRUE), "\n",
    "SVs soak up unwanted variation, but every SV costs a degree of freedom and an\n",
    "SV can also absorb real biology if it is correlated with your factor.\n",
    "Always check what the SVs correlate with before trusting them.\n")

# removeBatchEffect: ONLY for visualization
vsd <- vst(dds, blind = FALSE)
mat_clean <- removeBatchEffect(assay(vsd), covariates = sv,
                               design = model.matrix(~ sex + time, cd))
pc_raw <- prcomp(t(assay(vsd)[top, ])); pc_cl <- prcomp(t(mat_clean[top, ]))
pp <- rbind(data.frame(pc_raw$x[, 1:2], cd[, c("time", "sex")], which = "VST"),
            data.frame(pc_cl$x[, 1:2], cd[, c("time", "sex")], which = "VST, SVs removed"))
p <- ggplot(pp, aes(PC1, PC2, colour = time, shape = sex)) + geom_point(size = 3) +
  facet_wrap(~ which, scales = "free") + labs(title = "removeBatchEffect -- for plots only")
save_plot(p, file.path(od, "pca_sv_removed.png"), width = 10, height = 4.5)
cat("Never run DE on removeBatchEffect output: the model would not know degrees\n",
    "of freedom were spent. Put batch/SVs in the DESIGN instead (as above).\n")
