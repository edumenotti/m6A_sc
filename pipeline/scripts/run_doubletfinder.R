#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args) || grepl("^--", args[[i + 1]])) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

get_df_fun <- function(name) {
  ns <- asNamespace("DoubletFinder")
  if (exists(name, envir = ns, inherits = FALSE)) {
    return(get(name, envir = ns))
  }
  v3_name <- paste0(name, "_v3")
  if (exists(v3_name, envir = ns, inherits = FALSE)) {
    return(get(v3_name, envir = ns))
  }
  stop("DoubletFinder function not found: ", name, " or ", v3_name)
}

run_doubletfinder_compat <- function(seu, PCs, pN, pK, nExp) {
  real.cells <- rownames(seu@meta.data)
  counts <- if (SeuratObject::Version(seu) >= "5.0") {
    LayerData(seu, assay = "RNA", layer = "counts")
  } else {
    GetAssayData(object = seu, assay = "RNA", slot = "counts")
  }
  data <- counts[, real.cells]
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells / (1 - pN) - n_real.cells)
  message("Creating ", n_doublets, " artificial doublets with DoubletFinder-compatible scoring.")
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
  colnames(doublets) <- paste0("X", seq_len(n_doublets))
  data_wdoublets <- cbind(data, doublets)

  seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
  seu_wdoublets <- NormalizeData(seu_wdoublets, verbose = FALSE)
  seu_wdoublets <- FindVariableFeatures(
    seu_wdoublets,
    selection.method = "vst",
    nfeatures = min(2000, nrow(seu_wdoublets)),
    verbose = FALSE
  )
  seu_wdoublets <- ScaleData(seu_wdoublets, features = VariableFeatures(seu_wdoublets), verbose = FALSE)
  seu_wdoublets <- RunPCA(seu_wdoublets, features = VariableFeatures(seu_wdoublets), npcs = max(PCs), verbose = FALSE)
  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs, drop = FALSE]
  cell.names <- rownames(seu_wdoublets@meta.data)
  nCells <- length(cell.names)
  rm(seu_wdoublets)
  gc()

  dist.mat <- fields::rdist(pca.coord)
  pANN <- numeric(n_real.cells)
  names(pANN) <- real.cells
  k <- max(1, round(nCells * pK))
  k <- min(k, nCells - 1)
  for (i in seq_len(n_real.cells)) {
    neighbors <- order(dist.mat[, i])
    neighbors <- neighbors[2:(k + 1)]
    pANN[[i]] <- sum(neighbors > n_real.cells) / k
  }

  classifications <- rep("Singlet", n_real.cells)
  classifications[order(pANN, decreasing = TRUE)[seq_len(min(nExp, n_real.cells))]] <- "Doublet"
  seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data)]
  seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
  seu
}

install_doubletfinder <- function() {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      stop("DoubletFinder is missing and r-remotes/remotes is not installed.")
    }
    remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = "never")
  }
}

args <- parse_args()
if (isTRUE(args[["install-missing"]])) {
  install_doubletfinder()
}
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  stop("DoubletFinder is not installed. Re-run with --install-missing or install it with remotes::install_github().")
}

counts <- Matrix::readMM(args$counts)
genes <- readLines(args$genes, warn = FALSE)
cells <- readLines(args$cells, warn = FALSE)
rownames(counts) <- make.unique(genes)
colnames(counts) <- cells

expected_rate <- as.numeric(args$expected_rate)
pcs <- seq_len(as.integer(args$pcs))
seed <- as.integer(args$seed)
set.seed(seed)

obj <- CreateSeuratObject(counts = counts, project = args$sample_id, min.cells = 0, min.features = 0)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(2000, nrow(obj)), verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = max(pcs), verbose = FALSE)
obj <- FindNeighbors(obj, dims = pcs, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

paramSweep <- get_df_fun("paramSweep")
summarizeSweep <- get_df_fun("summarizeSweep")
find.pK <- get_df_fun("find.pK")
modelHomotypic <- get_df_fun("modelHomotypic")

auto_pk <- isTRUE(args[["auto-pk"]])
if (auto_pk) {
  sweep_res <- paramSweep(obj, PCs = pcs, sct = FALSE, num.cores = as.integer(args$num_cores))
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  pk_table <- find.pK(sweep_stats)
  pk_table$pK_numeric <- as.numeric(as.character(pk_table$pK))
  best <- pk_table[which.max(pk_table$BCmetric), , drop = FALSE]
  pK <- best$pK_numeric[[1]]
  if (is.na(pK) || pK <= 0) {
    pK <- as.numeric(args$pk)
  }
} else {
  pK <- as.numeric(args$pk)
}

n_exp <- max(1, round(expected_rate * ncol(obj)))
homotypic_prop <- modelHomotypic(obj$seurat_clusters)
n_exp_adj <- max(1, round(n_exp * (1 - homotypic_prop)))

obj <- run_doubletfinder_compat(
  obj,
  PCs = pcs,
  pN = 0.25,
  pK = pK,
  nExp = n_exp_adj
)

meta <- obj@meta.data
class_col <- grep("^DF.classifications", colnames(meta), value = TRUE)
pann_col <- grep("^pANN", colnames(meta), value = TRUE)
if (length(class_col) == 0 || length(pann_col) == 0) {
  stop("DoubletFinder did not return expected pANN/classification columns.")
}

calls <- data.frame(
  barcode = rownames(meta),
  doubletfinder_score = meta[[tail(pann_col, 1)]],
  doubletfinder_class = meta[[tail(class_col, 1)]],
  pK = pK,
  n_exp = n_exp,
  n_exp_adjusted = n_exp_adj,
  homotypic_prop = homotypic_prop,
  stringsAsFactors = FALSE
)
write.csv(calls, args$out_csv, row.names = FALSE)

summary <- data.frame(
  sample_id = args$sample_id,
  n_cells = ncol(obj),
  pK = pK,
  n_exp = n_exp,
  n_exp_adjusted = n_exp_adj,
  homotypic_prop = homotypic_prop,
  n_doublets = sum(calls$doubletfinder_class == "Doublet"),
  stringsAsFactors = FALSE
)
write.csv(summary, args$summary_csv, row.names = FALSE)
