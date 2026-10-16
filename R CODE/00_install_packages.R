# =============================================================================
# 00_install_packages.R  --  one-time installation of everything the course uses
# =============================================================================
# Needs R >= 4.3. Run once:   Rscript 00_install_packages.R
# Bioconductor packages are installed with BiocManager so that versions match
# your R release. Re-running is harmless: installed packages are skipped.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "dplyr", "tidyr",
          "matrixStats", "patchwork", "data.table", "WGCNA")
bioc <- c("DESeq2", "edgeR", "limma", "apeglm", "sva", "SummarizedExperiment",
          "org.Hs.eg.db", "GO.db", "AnnotationDbi", "fgsea", "BiocParallel",
          "genefilter", "airway")
optional <- c("msigdbr", "clusterProfiler", "org.Mm.eg.db", "IHW", "ashr")

need <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]

todo <- need(c(cran, bioc))
if (length(todo)) BiocManager::install(todo, update = FALSE, ask = FALSE)

# Optional extras used in "going further" sections; failures are not fatal
opt <- need(optional)
if (length(opt)) try(BiocManager::install(opt, update = FALSE, ask = FALSE))

missing <- need(c(cran, bioc))
if (length(missing)) {
  stop("These packages failed to install: ", paste(missing, collapse = ", "))
} else {
  message("\nAll core packages installed. R ", getRversion(),
          ", Bioconductor ", BiocManager::version())
}
