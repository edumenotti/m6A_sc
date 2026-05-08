#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}
has_flag <- function(flag) flag %in% args
if ("--help" %in% args || "-h" %in% args) {
  cat(paste(
    "Usage: prepare_nestorowa_reference.R [--out PATH] [--label-col COL|auto]",
    "[--method zellkonverter|mtx|both] [--anndata-version VERSION] [--install-missing]\n"
  ))
  quit(status = 0)
}
opt <- list(
  out = get_arg("--out", "data/references/nestorowa_2016_hspc.h5ad"),
  label_col = get_arg("--label-col", "auto"),
  method = get_arg("--method", "zellkonverter"),
  anndata_version = get_arg("--anndata-version", "0.10.2"),
  install_missing = has_flag("--install-missing")
)
if (!opt$method %in% c("zellkonverter", "mtx", "both")) {
  stop("--method must be one of: zellkonverter, mtx, both")
}

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0L) {
    return(invisible(NULL))
  }
  if (!opt$install_missing) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Re-run with --install-missing or install them in the pixi environment.")
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}

install_if_missing(c(
  "scRNAseq",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "Matrix",
  "AnnotationDbi",
  "org.Mm.eg.db",
  if (opt$method %in% c("zellkonverter", "both")) "zellkonverter" else character(0)
))

suppressPackageStartupMessages({
  library(scRNAseq)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(Matrix)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  if (opt$method %in% c("zellkonverter", "both")) {
    library(zellkonverter)
  }
})

map_ensembl_to_symbols <- function(sce) {
  gene_ids <- rownames(sce)
  if (is.null(gene_ids)) {
    stop("Nestorowa reference has no rownames/gene identifiers.")
  }
  if (!any(grepl("^ENSMUSG", gene_ids))) {
    rowData(sce)$gene_id <- gene_ids
    rowData(sce)$gene_symbol <- gene_ids
    return(sce)
  }

  symbols <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = gene_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  symbols <- as.character(symbols)
  keep <- !is.na(symbols) & symbols != "" & !duplicated(symbols)
  message("Mapped Ensembl IDs to gene symbols: keeping ", sum(keep), " of ", length(gene_ids),
          " genes; dropped ", sum(is.na(symbols) | symbols == ""), " unmapped and ",
          sum(duplicated(symbols[!is.na(symbols) & symbols != ""])), " duplicate symbols.")

  sce <- sce[keep, ]
  rowData(sce)$gene_id <- gene_ids[keep]
  rowData(sce)$gene_symbol <- symbols[keep]
  rownames(sce) <- symbols[keep]
  sce
}

message("Loading NestorowaHSCData from scRNAseq...")
sce <- scRNAseq::NestorowaHSCData(remove.htseq = TRUE, location = FALSE)
sce <- map_ensembl_to_symbols(sce)

col_df <- as.data.frame(S4Vectors::DataFrame(colData(sce)))
candidate_labels <- c(
  "cell_type", "cell.type", "CellType", "celltype", "cell_type1",
  "cell_type2", "label", "Label", "cell_label", "Cell.type", "mapping"
)

label_col <- opt$label_col
if (identical(label_col, "auto")) {
  present <- candidate_labels[candidate_labels %in% colnames(col_df)]
  if (length(present) == 0L) {
    message("Available colData columns:")
    message(paste(colnames(col_df), collapse = ", "))
    stop("Could not auto-detect a cell type label column. Re-run with --label-col.")
  }
  label_col <- present[[1]]
}
if (!label_col %in% colnames(col_df)) {
  stop("Label column not found: ", label_col)
}

labels <- as.character(col_df[[label_col]])
keep <- !is.na(labels) & labels != "" & labels != "NA"
sce <- sce[, keep]
labels <- labels[keep]

original_cells <- colnames(sce)
if (is.null(original_cells)) {
  original_cells <- paste0("cell_", seq_len(ncol(sce)))
}

colData(sce) <- S4Vectors::DataFrame(
  cell_type = labels,
  reference_batch = "Nestorowa2016",
  original_cell = original_cells,
  reference = "Nestorowa2016"
)
colnames(sce) <- paste0("Nestorowa2016_", make.unique(original_cells))

if (opt$method %in% c("zellkonverter", "both")) {
  out_dir <- dirname(opt$out)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  message("Writing H5AD with zellkonverter: ", opt$out)
  zellkonverter::writeH5AD(
    sce,
    opt$out,
    X_name = "counts",
    version = opt$anndata_version
  )
}

if (!opt$method %in% c("mtx", "both")) {
  message("Done.")
  message("Cells: ", ncol(sce), "; genes: ", nrow(sce), "; labels: ", length(unique(labels)))
  quit(status = 0)
}

export_dir <- paste0(tools::file_path_sans_ext(opt$out), "_export")
if (!dir.exists(export_dir)) {
  dir.create(export_dir, recursive = TRUE)
}

counts <- SummarizedExperiment::assay(sce, "counts")
if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

message("Writing Matrix Market export: ", export_dir)
Matrix::writeMM(counts, file.path(export_dir, "counts.genes_by_cells.mtx"))
write.table(
  data.frame(gene_id = rownames(sce), gene_symbol = rownames(sce)),
  file = file.path(export_dir, "genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.csv(
  as.data.frame(colData(sce)),
  file = file.path(export_dir, "obs.csv"),
  quote = TRUE,
  row.names = TRUE
)
writeLines(opt$out, file.path(export_dir, "target_h5ad.txt"))
message("Done.")
message("Cells: ", ncol(sce), "; genes: ", nrow(sce), "; labels: ", length(unique(labels)))
if (identical(opt$method, "mtx")) {
  message("Next: pixi run python pipeline/scripts/reference_export_to_h5ad.py --export-dir ",
          export_dir, " --out ", opt$out)
}
