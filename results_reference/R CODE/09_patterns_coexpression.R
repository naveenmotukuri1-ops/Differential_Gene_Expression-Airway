# =============================================================================
# MODULE 09  --  Expression patterns and co-expression networks
#                (GSE96870 mouse time course; needs module 07 results)
# =============================================================================
# Learning goals
#   * Turning thousands of LRT-significant genes into a few temporal patterns
#   * Clustering on scaled VST values; choosing k; reading cluster profiles
#   * WGCNA: soft thresholding, modules, eigengenes, module-trait correlation,
#     hub genes -- and when NOT to use WGCNA
#
# Run:  Rscript 09_patterns_coexpression.R   Outputs: results/09_patterns/
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
})
od <- out_dir("09_patterns")
set.seed(7)

f07 <- file.path(RESULTS_DIR, "07_multifactor", "mouse_dds.rds")
if (!file.exists(f07)) stop("Run 07_multifactor_designs.R first")
m7 <- readRDS(f07)
dds <- m7$dds; res_lrt <- m7$res_lrt; cd <- m7$coldata

# -----------------------------------------------------------------------------
banner("1. Pick the genes that change over time")
# -----------------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
sig <- rownames(res_lrt)[which(res_lrt$padj < 0.01)]
# Keep the strongest 2000 to make clustering readable
sig <- head(sig[order(res_lrt[sig, "padj"])], 2000)
cat("Genes used for pattern discovery:", length(sig), "\n")

# Remove the sex effect from the VALUES used for clustering (visual/exploratory
# use of removeBatchEffect is fine; testing was already done with sex in the model)
mat <- limma::removeBatchEffect(assay(vsd)[sig, ], batch = cd$sex,
                                design = model.matrix(~ time, cd))
z <- t(scale(t(mat)))               # z-score per gene: we cluster SHAPES, not levels

# -----------------------------------------------------------------------------
banner("2. How many patterns? (k-means, within-cluster SS)")
# -----------------------------------------------------------------------------
# Average per time point first: patterns are about time, and averaging cuts noise
zt <- sapply(levels(cd$time), function(t) rowMeans(z[, cd$time == t, drop = FALSE]))
wss <- sapply(2:8, function(k) kmeans(zt, k, nstart = 25)$tot.withinss)
save_plot(function() plot(2:8, wss, type = "b", xlab = "k", ylab = "total within SS",
                          main = "Elbow plot"), file.path(od, "elbow.png"), 5, 4)
k <- 4
cat("Using k =", k, "(with 3 time points there are only a few possible shapes:\n",
    " up-early, up-late, down-early, down-late, transient peak/dip)\n")
km <- kmeans(zt, k, nstart = 50)
cl <- km$cluster
print(table(cluster = cl))

prof <- data.frame(gene = rep(rownames(zt), ncol(zt)),
                   time = factor(rep(colnames(zt), each = nrow(zt)), levels = levels(cd$time)),
                   z = c(zt), cluster = paste0("C", rep(cl, ncol(zt))))
prof$cluster <- paste0(prof$cluster, " (n=", table(cl)[sub("C", "", prof$cluster)], ")")
p <- ggplot(prof, aes(time, z, group = gene)) + geom_line(alpha = .05) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", colour = "red", linewidth = 1.2) +
  facet_wrap(~ cluster) + labs(y = "z-scored expression", title = "Temporal patterns")
save_plot(p, file.path(od, "cluster_profiles.png"), width = 8, height = 6)

ord <- order(cl)
ann <- data.frame(time = cd$time, sex = cd$sex, row.names = colnames(z))
rann <- data.frame(cluster = factor(paste0("C", cl)), row.names = rownames(z))
cols <- order(cd$time, cd$sex)
save_plot(function() pheatmap(z[ord, cols], cluster_rows = FALSE, cluster_cols = FALSE,
                              annotation_col = ann, annotation_row = rann, show_rownames = FALSE,
                              breaks = seq(-2.5, 2.5, length.out = 101),
                              main = "LRT genes by temporal cluster"),
          file.path(od, "cluster_heatmap.png"), width = 8, height = 8)

# Describe each cluster: its centre and its most significant genes
for (i in seq_len(k)) {
  g <- names(cl)[cl == i]
  g <- g[order(res_lrt[g, "padj"])]
  cat("C", i, " centre (Day0,Day4,Day8) = ", paste(round(km$centers[i, ], 2), collapse = ", "),
      "\n   most significant: ", paste(head(g, 10), collapse = ", "), "\n", sep = "")
}
cat("Two clusters may look alike (e.g. 'late up' and 'gradual up'): k-means always\n",
    "returns k groups. Try k = 3, 5, 6 and keep the k whose shapes you can name.\n")
write.csv(data.frame(gene = names(cl), cluster = cl, padj_LRT = res_lrt[names(cl), "padj"]),
          file.path(od, "gene_clusters.csv"), row.names = FALSE)
cat("Next step (exercise): run ORA per cluster with org.Mm.eg.db as in module 06,\n",
    "using all LRT-tested genes as the universe.\n")

# -----------------------------------------------------------------------------
banner("3. WGCNA co-expression modules")
# -----------------------------------------------------------------------------
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  cat("WGCNA not installed -- skipping this section.\n")
} else {
  suppressPackageStartupMessages(library(WGCNA))
  allowWGCNAThreads(2)
  # WGCNA wants samples in rows; use the most variable genes (not only DE genes --
  # selecting on DE then 'discovering' trait-associated modules is circular)
  vv <- assay(vsd)
  top <- head(order(apply(vv, 1, var), decreasing = TRUE), 3000)
  datExpr <- t(vv[top, ])
  stopifnot(goodSamplesGenes(datExpr, verbose = 0)$allOK)

  sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, 12, 14, 16, 18, 20),
                           networkType = "signed", verbose = 0)
  fi <- sft$fitIndices
  fit_r2 <- -sign(fi$slope) * fi$SFT.R.sq
  power <- fi$Power[which(fit_r2 > 0.8)[1]]
  if (is.na(power)) {
    # The WGCNA FAQ: if no power reaches ~0.8, a strong global driver (here, the
    # infection time course) usually dominates. Its recommended fallback for a
    # SIGNED network with 20-30 samples is power 16.
    power <- 16
    cat("Scale-free fit never reached 0.8 -- typical when one factor (time) drives\n",
        "thousands of genes. Using the WGCNA FAQ fallback for signed networks.\n")
  }
  cat("Chosen soft-threshold power:", power, "\n")
  save_plot(function() {
    plot(fi$Power, fit_r2, type = "n", xlab = "soft threshold power",
         ylab = "scale-free fit (signed R^2)", main = "Pick the lowest power reaching ~0.8")
    text(fi$Power, fit_r2, fi$Power, col = "red"); abline(h = .8, col = "blue")
  }, file.path(od, "wgcna_soft_threshold.png"), 6, 5)

  net <- blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed",
                          minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = FALSE,
                          maxBlockSize = 5000, verbose = 0)
  cat("Modules (grey = unassigned):\n"); print(table(net$colors))

  traits <- data.frame(Day4 = as.numeric(cd$time == "Day4"),
                       Day8 = as.numeric(cd$time == "Day8"),
                       infected = as.numeric(cd$infection == "InfluenzaA"),
                       male = as.numeric(cd$sex == "Male"))
  MEs <- orderMEs(net$MEs)
  mt_cor <- cor(MEs, traits, use = "p")
  mt_p <- corPvalueStudent(mt_cor, nrow(datExpr))
  txt <- matrix(paste0(round(mt_cor, 2), "\n(", signif(mt_p, 1), ")"), nrow(mt_cor))
  save_plot(function() {
    par(mar = c(6, 9, 3, 3))
    labeledHeatmap(mt_cor, xLabels = colnames(traits), yLabels = rownames(mt_cor),
                   ySymbols = rownames(mt_cor), colorLabels = FALSE,
                   colors = blueWhiteRed(50), textMatrix = txt, setStdMargins = FALSE,
                   cex.text = 0.6, zlim = c(-1, 1), main = "Module eigengene - trait correlation")
  }, file.path(od, "wgcna_module_trait.png"), 7, 7)
  print(round(mt_cor, 2))

  # Hub genes: highest module membership (kME) in the module most tied to Day8
  best <- rownames(mt_cor)[which.max(abs(mt_cor[, "Day8"]))]
  mod <- sub("^ME", "", best)
  kME <- cor(datExpr[, net$colors == mod], MEs[, best])
  hubs <- head(sort(kME[, 1], decreasing = TRUE), 15)
  cat("\nModule most associated with Day8:", mod, " (r =", round(mt_cor[best, "Day8"], 2), ")\n",
      "Hub genes (highest kME):", paste(names(hubs), collapse = ", "), "\n")
  # Where did the biology from module 07 land?
  markers <- c("Irf7", "Isg15", "Ifit1", "Ifit3", "Oasl2", "Bst2", "Plp1", "Mog")
  markers <- markers[markers %in% colnames(datExpr)]
  cat("\nModule of known marker genes:\n")
  print(setNames(net$colors[match(markers, colnames(datExpr))], markers))
  isg_mod <- net$colors[match("Irf7", colnames(datExpr))]
  if (!is.na(isg_mod) && isg_mod != "grey") {
    cat("The interferon genes share one module ('", isg_mod, "', Day4 r = ",
        round(mt_cor[paste0("ME", isg_mod), "Day4"], 2), "): a TRANSIENT program that\n",
        "the k-means step above did not isolate, because the LRT top-2000 list is\n",
        "dominated by the larger Day8 changes. Different lenses find different things.\n", sep = "")
  }
  write.csv(data.frame(gene = colnames(datExpr), module = net$colors),
            file.path(od, "wgcna_modules.csv"), row.names = FALSE)
  cat("\nWGCNA caveats: needs >= ~15-20 samples (we have 22), is sensitive to\n",
      "outliers and batch, and module-trait correlations are exploratory. A hub gene\n",
      "is a hypothesis for follow-up, not a validated regulator.\n")
}
