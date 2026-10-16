# =============================================================================
# MODULE 06  --  From gene lists to biology: ORA, GSEA, camera   (airway)
# =============================================================================
# Learning goals
#   * Over-representation analysis (ORA) = a hypergeometric test, done by hand
#   * The background (universe) mistake and how much it inflates results
#   * GSEA with fgsea on a ranked list (no arbitrary cutoff)
#   * Competitive gene-set testing with limma::camera (accounts for gene-gene
#     correlation -- ORA and simple GSEA do not)
#   * Collapsing redundant GO terms; reading results critically
#
# Gene sets come from GO via org.Hs.eg.db + GO.db, so this module works offline.
# If msigdbr is installed, a Hallmark section runs too.
#
# Run:  Rscript 06_functional_enrichment.R    Outputs: results/06_enrichment/
# Needs 04 and 05 results.
# =============================================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(fgsea)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
})
od <- out_dir("06_enrichment")
set.seed(42)

de2 <- readRDS(file.path(RESULTS_DIR, "04_deseq2", "airway_deseq2.rds"))
lim <- readRDS(file.path(RESULTS_DIR, "05_compare", "airway_limma.rds"))
res <- de2$res

# -----------------------------------------------------------------------------
banner("1. Build GO Biological Process gene sets (Entrez IDs)")
# -----------------------------------------------------------------------------
tested <- rownames(res)[!is.na(res$padj)]
entrez <- ens_to_entrez(tested)
cat("Tested genes:", length(tested), " with Entrez ID:", sum(!is.na(entrez)), "\n")
# Keep one Ensembl gene per Entrez ID (highest baseMean) to avoid double counting
map <- data.frame(ens = tested, entrez = entrez, bm = res[tested, "baseMean"])
map <- map[!is.na(map$entrez), ]
map <- map[order(-map$bm), ]
map <- map[!duplicated(map$entrez), ]
universe <- map$entrez

go2genes <- as.list(org.Hs.egGO2ALLEGS)                 # GO term -> genes (incl. children)
ont <- suppressMessages(AnnotationDbi::select(GO.db, keys = names(go2genes),
                                              columns = c("ONTOLOGY", "TERM"), keytype = "GOID"))
bp_ids <- ont$GOID[ont$ONTOLOGY == "BP"]
go_bp <- lapply(go2genes[bp_ids], function(g) intersect(unique(g), universe))
sizes <- lengths(go_bp)
go_bp <- go_bp[sizes >= 15 & sizes <= 500]              # too small = noisy, too big = vague
term_name <- setNames(ont$TERM, ont$GOID)
cat("GO BP sets with 15-500 genes in our universe:", length(go_bp), "\n")

# -----------------------------------------------------------------------------
banner("2. ORA by hand: hypergeometric test")
# -----------------------------------------------------------------------------
up <- map$entrez[map$ens %in% rownames(res)[which(res$padj < .05 & res$log2FoldChange > 1)]]
cat("Up-regulated genes (padj<.05, LFC>1) with Entrez:", length(up), "\n")

ora <- function(hits, universe, sets) {
  N <- length(universe); n <- length(hits)
  out <- vapply(sets, function(s) {
    K <- length(s); k <- sum(hits %in% s)
    c(k = k, K = K, p = phyper(k - 1, K, N - K, n, lower.tail = FALSE))
  }, numeric(3))
  out <- as.data.frame(t(out))
  out$expected <- out$K * n / N
  out$fold <- out$k / out$expected
  out$padj <- p.adjust(out$p, "BH")
  out$term <- term_name[rownames(out)]
  out[order(out$p), c("term", "k", "K", "expected", "fold", "p", "padj")]
}
ora_up <- ora(up, universe, go_bp)
cat("Top terms for up-regulated genes:\n")
print(head(ora_up[, c("term", "k", "K", "fold", "padj")], 12), digits = 3)
write.csv(ora_up, file.path(od, "ORA_GO_BP_up.csv"))

# The background mistake: use ALL annotated human genes as universe
all_genes <- unique(unlist(go2genes[bp_ids]))
go_bp_all <- lapply(go2genes[names(go_bp)], unique)
ora_wrong <- ora(up, all_genes, go_bp_all)
ora_wrong <- ora_wrong[rownames(ora_up), ]
sig_ok <- ora_up$padj < .05; sig_bad <- ora_wrong$padj < .05
cat("\nSignificant terms (FDR<5%), correct universe (tested genes):", sum(sig_ok), "\n",
    "Significant terms (FDR<5%), WRONG universe (whole genome):  ", sum(sig_bad),
    "  (shared:", sum(sig_ok & sig_bad), ")\n")
shift <- data.frame(term = ora_up$term,
                    log10p_correct = round(log10(ora_up$p), 1),
                    log10p_wrong   = round(log10(ora_wrong$p), 1),
                    K_expressed = ora_up$K, K_genome = ora_wrong$K)
shift$ratio_K <- round(shift$K_genome / shift$K_expressed, 2)
shift <- shift[order(shift$log10p_wrong - shift$log10p_correct), ]
cat("\nTerms made MORE significant by the wrong universe (mostly-expressed sets):\n")
print(head(shift[shift$log10p_correct < -2 | shift$log10p_wrong < -2, ], 6), row.names = FALSE)
cat("\nTerms made LESS significant (sets full of genes not expressed here):\n")
print(head(shift[nrow(shift):1, ][shift$log10p_correct[nrow(shift):1] < -3, ], 6), row.names = FALSE)
cat("\nThe genome background changes the answer in BOTH directions, depending on\n",
    "how many of a term's genes are expressed in this tissue (ratio_K). The\n",
    "question you want answered is 'among genes that COULD have been called DE,\n",
    "which functions are over-represented?' -- so the universe is the tested genes.\n")
write.csv(shift, file.path(od, "ORA_universe_effect.csv"), row.names = FALSE)

# limma::goana does the same ORA (with the right universe) in one line
ga <- goana(list(Up = up), universe = universe, species = "Hs")
ga <- ga[ga$Ont == "BP", ]
cat("\nlimma::goana top 5 (should look like the hand-made ORA):\n")
print(head(topGO(ga, ontology = "BP", number = 5)), digits = 3)

# -----------------------------------------------------------------------------
banner("3. GSEA with fgsea (uses ALL genes, no cutoff)")
# -----------------------------------------------------------------------------
# Rank statistic: Wald stat combines effect size and confidence. (Shrunken LFC
# is another good choice; raw p-values are not -- they lose the direction.)
rank_stat <- setNames(res[map$ens, "stat"], map$entrez)
rank_stat <- sort(rank_stat[!is.na(rank_stat)], decreasing = TRUE)
gsea <- fgsea(pathways = go_bp, stats = rank_stat, minSize = 15, maxSize = 500, eps = 0)
gsea$term <- term_name[gsea$pathway]
gsea <- gsea[order(gsea$pval), ]
cat("GSEA significant (padj<5%):", sum(gsea$padj < .05),
    " (up:", sum(gsea$padj < .05 & gsea$NES > 0), " down:", sum(gsea$padj < .05 & gsea$NES < 0), ")\n")
print(head(gsea[, c("term", "NES", "size", "padj")], 15), digits = 3)

# Collapse redundant GO terms (parent/child terms share genes)
sig_paths <- gsea[gsea$padj < .01, ]
if (nrow(sig_paths) > 1) {
  coll <- collapsePathways(sig_paths, go_bp, rank_stat)
  main <- sig_paths[sig_paths$pathway %in% coll$mainPathways, ]
  cat("\nAfter collapsing redundant terms:", nrow(main), "of", nrow(sig_paths), "remain\n")
  print(head(main[order(main$pval), c("term", "NES", "padj")], 15), digits = 3)
}
gsea_out <- gsea
gsea_out$leadingEdge <- vapply(gsea_out$leadingEdge, paste, character(1), collapse = ";")
write.csv(gsea_out, file.path(od, "GSEA_GO_BP.csv"), row.names = FALSE)

top_up <- gsea$pathway[gsea$NES > 0][1]
save_plot(plotEnrichment(go_bp[[top_up]], rank_stat) +
            labs(title = term_name[top_up]),
          file.path(od, "gsea_top_up.png"), width = 7, height = 4)
topP <- c(head(gsea$pathway[gsea$NES > 0], 8), head(gsea$pathway[gsea$NES < 0], 8))
bar <- gsea[match(topP, gsea$pathway), ]
p <- ggplot(bar, aes(NES, reorder(substr(term, 1, 55), NES), fill = NES > 0)) + geom_col() +
  labs(y = NULL, title = "Top GO BP terms by GSEA") + theme(legend.position = "none")
save_plot(p, file.path(od, "gsea_top_terms.png"), width = 9, height = 6)

# -----------------------------------------------------------------------------
banner("4. camera: competitive test that respects gene-gene correlation")
# -----------------------------------------------------------------------------
v <- lim$v
v_entrez <- ens_to_entrez(rownames(v))
idx <- ids2indices(go_bp, v_entrez)
ctr <- which(colnames(lim$design) == "dextreated")
# Default: a small fixed inter-gene correlation (0.01) -- the recommended setting
cam <- camera(v, idx, design = lim$design, contrast = ctr)
cam$term <- term_name[rownames(cam)]
cat("camera (inter.gene.cor = 0.01) significant (FDR<5%):", sum(cam$FDR < .05), "\n")
print(head(cam[, c("term", "NGenes", "Direction", "FDR")], 10), digits = 3)
# Estimating the correlation from only 8 samples is very noisy -> very conservative
cam_est <- camera(v, idx, design = lim$design, contrast = ctr, inter.gene.cor = NULL)
cat("\ncamera with correlation ESTIMATED per set: significant =", sum(cam_est$FDR < .05), "\n",
    "Compare with fgsea's", sum(gsea$padj < .05), "terms. Genes in a pathway are\n",
    "co-regulated, so they are not independent evidence; tests that assume\n",
    "independence are over-confident. With 8 samples the correlation estimate is\n",
    "too noisy, hence the fixed 0.01 default. Treat GO results as hypotheses.\n")
write.csv(cam, file.path(od, "camera_GO_BP.csv"))

# -----------------------------------------------------------------------------
banner("5. (Optional) MSigDB Hallmark sets")
# -----------------------------------------------------------------------------
if (requireNamespace("msigdbr", quietly = TRUE)) {
  hm <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  hallmark <- split(as.character(hm$ncbi_gene), hm$gs_name)
  gh <- fgsea(hallmark, rank_stat, minSize = 10, maxSize = 500, eps = 0)
  gh <- gh[order(gh$pval), ]
  print(head(gh[, c("pathway", "NES", "padj")], 15), digits = 3)
  gh$leadingEdge <- vapply(gh$leadingEdge, paste, character(1), collapse = ";")
  write.csv(gh, file.path(od, "GSEA_hallmark.csv"), row.names = FALSE)
} else {
  cat("msigdbr not installed -- skipping. install.packages('msigdbr') to enable.\n")
}
