# =============================================================================
# 01_get_data.R  --  download / rebuild the two course datasets into data/
# =============================================================================
# The course ships with data/ already filled in, so you only need this script
# if data/ is missing or you want to see where the data comes from.
#
# Dataset 1  "airway"   (human, modules 02-06, 08)
#   Himes et al. 2014 PLoS One, GEO GSE52778. Four primary human airway
#   smooth-muscle cell lines (donors), each untreated or treated with
#   dexamethasone. 8 samples, paired design. Counts are the Bioconductor
#   `airway` package (STAR alignment + GenomicAlignments::summarizeOverlaps,
#   Ensembl GRCh37.75 genes).
#
# Dataset 2  "GSE96870" (mouse, modules 07, 09)
#   Blackmore et al. 2017 PNAS. Cerebellum of mice infected intranasally
#   with influenza A; Day0 (non-infected), Day4, Day8; both sexes. 22 samples.
#   Counts as curated for the Carpentries / Bioconductor RNA-seq lesson.
#
# Run from the course root:   Rscript 01_get_data.R
# =============================================================================

dir.create("data", showWarnings = FALSE)

# ---------------------------------------------------------------- airway -----
build_airway <- function() {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  if (requireNamespace("airway", quietly = TRUE)) {
    message("Using the installed Bioconductor 'airway' package")
    data("airway", package = "airway", envir = environment())
  } else {
    url <- "https://raw.githubusercontent.com/bioc/airway/devel/data/airway.RData"
    message("airway package not installed; downloading ", url)
    tmp <- tempfile(fileext = ".RData")
    download.file(url, tmp, mode = "wb", quiet = TRUE)
    load(tmp)
  }

  cts <- assay(airway, "counts")
  cd  <- as.data.frame(colData(airway))
  coldata <- data.frame(
    sample  = rownames(cd),
    cell    = as.character(cd$cell),                  # donor / cell line
    dex     = ifelse(cd$dex == "trt", "treated", "control"),
    geo     = as.character(cd$SampleName),
    avgLength = cd$avgLength,                         # mean read length
    stringsAsFactors = FALSE)

  # Gene annotation + exonic length (union of exons) -- needed for TPM in module 03
  rd <- as.data.frame(rowData(airway))
  exon_len <- sum(width(GenomicRanges::reduce(rowRanges(airway))))
  genes <- data.frame(gene_id = rd$gene_id, symbol = rd$gene_name,
                      biotype = rd$gene_biotype, chr = rd$seq_name,
                      exonic_length = as.integer(exon_len[rd$gene_id]),
                      stringsAsFactors = FALSE)

  out <- data.frame(gene_id = rownames(cts), cts, check.names = FALSE)
  con <- gzfile("data/airway_counts.csv.gz", "w")
  write.csv(out, con, row.names = FALSE); close(con)
  write.csv(coldata, "data/airway_coldata.csv", row.names = FALSE)
  con <- gzfile("data/airway_genes.csv.gz", "w")
  write.csv(genes, con, row.names = FALSE); close(con)
  message("airway: ", nrow(cts), " genes x ", ncol(cts), " samples written to data/")
}

# --------------------------------------------------------------- GSE96870 ----
build_mouse <- function() {
  base <- "https://raw.githubusercontent.com/Bioconductor/bioconductor-teaching/master/data/GSE96870/"
  files <- c("GSE96870_counts_cerebellum.csv", "GSE96870_coldata_cerebellum.csv",
             "GSE96870_rowranges.tsv")
  for (f in files) {
    dest <- file.path("data", f)
    if (file.exists(dest)) { message("  have ", dest); next }
    message("  downloading ", f)
    download.file(paste0(base, f), dest, mode = "wb", quiet = TRUE)
  }
  message("GSE96870: files in data/")
}

if (!file.exists("data/airway_counts.csv.gz")) build_airway() else message("airway data already present")
build_mouse()

# Sanity check: these MD5s are what the course was tested with
print(tools::md5sum(list.files("data", full.names = TRUE)))
