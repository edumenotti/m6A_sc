#!/usr/bin/env Rscript
# 21b_hemascribe.R
#
# Load the MTX bundle exported by 21a_export_for_hemascribe.py, build a
# Seurat object, run HemaScribe, and write per-cell labels (broad + fine)
# to CSV for 21c to merge back into the AnnData object.
#
# Usage:
#   pixi run Rscript pipeline/scripts/21b_hemascribe.R \
#     --indir  results/21_hemascribe/export \
#     --out    results/21_hemascribe/hemascribe_labels.csv

suppressPackageStartupMessages({
  library(argparse)
  library(Matrix)
  library(Seurat)
  library(HemaScribe)
})

p <- ArgumentParser()
p$add_argument("--indir", required = TRUE,
               help = "Directory with counts.mtx.gz / barcodes.tsv.gz / features.tsv.gz / obs.csv")
p$add_argument("--out", required = TRUE)
args <- p$parse_args()

dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)

mtx <- file.path(args$indir, "counts.mtx.gz")
bc  <- file.path(args$indir, "barcodes.tsv.gz")
ft  <- file.path(args$indir, "features.tsv.gz")
ob  <- file.path(args$indir, "obs.csv")
for (f in c(mtx, bc, ft, ob)) {
  if (!file.exists(f)) stop("missing: ", f)
}

cat("[21b] Reading MTX bundle ...\n")
counts <- ReadMtx(mtx = mtx, cells = bc, features = ft,
                  feature.column = 2, cell.column = 1)
cat("[21b] counts dim:", paste(dim(counts), collapse = " x "), "\n")

cat("[21b] Reading obs.csv ...\n")
meta <- read.csv(ob, row.names = 1, check.names = FALSE)
stopifnot(all(colnames(counts) %in% rownames(meta)))
meta <- meta[colnames(counts), , drop = FALSE]

cat("[21b] Building Seurat object ...\n")
seu <- CreateSeuratObject(counts = counts, meta.data = meta,
                          assay = "RNA", min.cells = 0, min.features = 0)
rm(counts); gc(verbose = FALSE)

cat("[21b] NormalizeData ...\n")
seu <- NormalizeData(seu, verbose = FALSE)

cat("[21b] Ensuring HemaScribe reference data is cached ...\n")
try(download_data(),        silent = TRUE)
try(download_mapping_data(), silent = TRUE)

cat("[21b] Running HemaScribe(return.full = FALSE) ...\n")
t0 <- Sys.time()
seu <- HemaScribe(seu, return.full = FALSE)
dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
cat("[21b] HemaScribe finished in", dt, "min\n")

md <- seu@meta.data
cat("[21b] HemaScribe-related meta columns:\n")
print(grep("hema|broad|fine|score|label",
           colnames(md), ignore.case = TRUE, value = TRUE))

pick <- function(md, patterns) {
  for (pat in patterns) {
    hit <- grep(pat, colnames(md), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}
broad_col <- pick(md, c("^broad\\.annot$", "^broad_annot$",
                        "^broad_classify$", "broad.?class", "broad.?label",
                        "broad.?annot"))
fine_col  <- pick(md, c("^fine\\.annot$", "^fine_annot$",
                        "^fine_classify$", "fine.?class", "fine.?label",
                        "fine.?annot"))
score_col <- pick(md, c("hematopoietic.?score", "hema.*score"))

cat(sprintf("[21b] resolved: broad=%s | fine=%s | score=%s\n",
            broad_col, fine_col, score_col))

if (is.na(broad_col) && is.na(fine_col)) {
  stop("HemaScribe did not produce broad/fine label columns. Inspect md.")
}

out <- data.frame(
  cell_id          = rownames(md),
  hemascribe_broad = if (!is.na(broad_col)) as.character(md[[broad_col]]) else NA_character_,
  hemascribe_fine  = if (!is.na(fine_col))  as.character(md[[fine_col]])  else NA_character_,
  hemascribe_score = if (!is.na(score_col)) as.numeric(md[[score_col]])    else NA_real_,
  stringsAsFactors = FALSE
)
write.csv(out, args$out, row.names = FALSE)
cat("[21b] Wrote", args$out, "with", nrow(out), "rows\n")

full_out <- sub("\\.csv$", "_full_meta.csv", args$out)
write.csv(cbind(cell_id = rownames(md), md), full_out, row.names = FALSE)
cat("[21b] Wrote", full_out, "(full meta)\n")
