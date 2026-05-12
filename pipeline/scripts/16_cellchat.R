#!/usr/bin/env Rscript
# 16_cellchat.R
#
# Cell-cell communication analysis using CellChat v2.
#
# For each condition (WT_DMSO, WT_STM, Mutant_DMSO, Mutant_STM) and for each
# of two annotation levels (manual_level1, manual_level2) builds one CellChat
# object, then performs cross-condition comparison (rankNet pathway info-flow,
# subsetCommunication for L-R differential, plus delta interaction counts).
#
# Input bundle (MatrixMarket export from export_h5ad_for_r.py)
#   counts.mtx, barcodes.txt, features.csv, metadata.csv, export_info.json
#
# Inputs
#   --input-dir  directory holding the MatrixMarket bundle
#   --out        output directory
#   --organism   "Mm" (default) or "Hs"
#   --min-cells  minimum cells per (celltype, condition) cell group to keep
#                (default 10 — matches filterCommunication default)
#
# Outputs (in --out/<level>/)
#   cellchat_<condition>.rds                  one CellChat object per condition
#   cellchat_interaction_counts.csv           #links per condition (total)
#   cellchat_interaction_weights.csv          total interaction strength per condition
#   cellchat_chord_<condition>.png            chord diagrams (counts; width-capped)
#   cellchat_pathway_infoflow.csv             rankNet output (info-flow per pathway, per condition)
#   cellchat_delta_interaction_counts.csv     pairwise sender×receiver delta
#   cellchat_diff_LR_<A_vs_B>.csv             differential L-R interactions per comparison
#   cellchat_merged.rds                       merged object across conditions
#   cellchat_celltype_counts.csv              cell counts per (celltype, condition)

suppressPackageStartupMessages({
  library(argparse)
  library(Matrix)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
})

# ── CLI ────────────────────────────────────────────────────────────────────
parser <- ArgumentParser()
parser$add_argument("--input-dir", required=TRUE,
                    help="dir with counts.mtx, barcodes.txt, features.csv, metadata.csv, export_info.json")
parser$add_argument("--out",       default=".")
parser$add_argument("--organism",  default="Mm", choices=c("Mm","Hs"))
parser$add_argument("--min-cells", type="integer", default=10)
args <- parser$parse_args()

dir.create(args$out, recursive=TRUE, showWarnings=FALSE)
options(future.globals.maxSize = 4 * 1024^3)  # 4 GB

# ── Load MatrixMarket bundle ───────────────────────────────────────────────
message("Loading expression bundle from ", args$input_dir, " ...")
info <- if (file.exists(file.path(args$input_dir, "export_info.json"))) {
  jsonlite::fromJSON(file.path(args$input_dir, "export_info.json"))
} else list(layer="unknown", integer=NA)
message("  layer='", info$layer, "', integer=", info$integer)

counts   <- readMM(file.path(args$input_dir, "counts.mtx"))
counts   <- as(counts, "CsparseMatrix")
barcodes <- readLines(file.path(args$input_dir, "barcodes.txt"))
features <- read.csv(file.path(args$input_dir, "features.csv"), row.names=1, check.names=FALSE)
meta     <- read.csv(file.path(args$input_dir, "metadata.csv"), row.names=1, check.names=FALSE)

rownames(counts) <- rownames(features)
colnames(counts) <- barcodes

# Hard guarantee that counts and meta are aligned cell-wise
stopifnot(ncol(counts) == nrow(meta))
stopifnot(identical(colnames(counts), rownames(meta)))

meta$condition <- paste(meta$genotype, meta$treatment, sep="_")

# CellChat expects log-normalized data. If we got raw integers, normalize now.
data_is_integer <- isTRUE(info$integer) ||
  (!is.na(info$integer) && info$integer == "true") ||
  all(counts@x == round(counts@x))
if (data_is_integer) {
  message("  Input looks like raw counts → applying CellChat::normalizeData()")
  data_input <- normalizeData(counts)
} else {
  message("  Input already normalized → skipping normalizeData()")
  data_input <- counts
}
rm(counts); gc(verbose=FALSE)

# ── CellChatDB ─────────────────────────────────────────────────────────────
CellChatDB.use <- if (args$organism == "Mm") CellChatDB.mouse else CellChatDB.human

# ── Helper: run one CellChat object (assumes data is normalized) ───────────
run_cellchat <- function(data_sub, meta_sub, label, min_cells) {
  message("  Building CellChat for: ", label,
          " (", ncol(data_sub), " cells, ",
          length(unique(meta_sub$celltype)), " cell types)")
  cc <- createCellChat(object=data_sub, meta=meta_sub, group.by="celltype")
  cc@DB <- CellChatDB.use
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type="triMean", population.size=TRUE)
  cc <- filterCommunication(cc, min.cells=min_cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# ── Per-level pipeline ─────────────────────────────────────────────────────
run_level <- function(level_col, out_subdir) {
  message("\n=== Annotation level: ", level_col, " ===")
  dir.create(out_subdir, recursive=TRUE, showWarnings=FALSE)

  meta_lvl <- meta
  meta_lvl$celltype <- as.character(meta_lvl[[level_col]])

  # Drop NA celltype cells
  keep_cell <- !is.na(meta_lvl$celltype) & meta_lvl$celltype != ""
  if (!all(keep_cell)) {
    message("  Dropping ", sum(!keep_cell), " cells with NA/empty ", level_col)
  }

  # Cell-count diagnostic table (celltype × condition) ─ before any filtering
  ct_table <- as.data.frame.matrix(table(
    celltype  = meta_lvl$celltype[keep_cell],
    condition = meta_lvl$condition[keep_cell]
  ))
  ct_table$celltype <- rownames(ct_table)
  ct_table <- ct_table[, c("celltype", setdiff(colnames(ct_table), "celltype"))]
  write.csv(ct_table, file.path(out_subdir, "cellchat_celltype_counts.csv"),
            row.names=FALSE)
  message("  Cell counts per (celltype × condition):")
  print(ct_table)

  conditions <- c("WT_DMSO", "WT_STM", "Mutant_DMSO", "Mutant_STM")
  conditions <- conditions[conditions %in% unique(meta_lvl$condition)]
  cc_list <- list()

  for (cond in conditions) {
    idx <- keep_cell & meta_lvl$condition == cond
    if (sum(idx) < 50) {
      message("  Skipping ", cond, " — only ", sum(idx), " cells")
      next
    }

    # Pre-filter small celltype groups BEFORE createCellChat
    ct_counts <- table(meta_lvl$celltype[idx])
    keep_ct <- names(ct_counts)[ct_counts >= args$min_cells]
    dropped <- setdiff(names(ct_counts), keep_ct)
    if (length(dropped) > 0) {
      message("  [", cond, "] dropping cell types with < ", args$min_cells,
              " cells: ", paste(dropped, collapse=", "))
    }
    idx <- idx & meta_lvl$celltype %in% keep_ct
    if (length(unique(meta_lvl$celltype[idx])) < 2) {
      message("  Skipping ", cond, " — fewer than 2 retained cell types")
      next
    }

    data_sub <- data_input[, idx]
    meta_sub <- meta_lvl[idx, , drop=FALSE]
    stopifnot(ncol(data_sub) == nrow(meta_sub))

    cc <- run_cellchat(data_sub, meta_sub, cond, args$min_cells)
    cc_list[[cond]] <- cc
    saveRDS(cc, file.path(out_subdir, paste0("cellchat_", cond, ".rds")))

    # Chord diagram with capped edge width so weighted lines stay legible
    png(file.path(out_subdir, paste0("cellchat_chord_", cond, ".png")),
        width=1200, height=1200, res=150)
    netVisual_circle(cc@net$count,
                     weight.scale = TRUE,
                     edge.weight.max = max(cc@net$count, na.rm=TRUE),
                     edge.width.max = 8,
                     label.edge = FALSE,
                     title.name = paste0("Interaction counts — ", cond))
    dev.off()
  }

  if (length(cc_list) < 2) {
    message("  Fewer than 2 conditions retained — skipping cross-condition section")
    return(invisible(NULL))
  }

  # ── Cross-condition summary tables ──────────────────────────────────────
  summary_df <- do.call(rbind, lapply(names(cc_list), function(nm) {
    data.frame(
      condition    = nm,
      total_links  = sum(cc_list[[nm]]@net$count > 0),
      total_weight = sum(cc_list[[nm]]@net$weight, na.rm=TRUE)
    )
  }))
  write.csv(summary_df[, c("condition","total_links")],
            file.path(out_subdir, "cellchat_interaction_counts.csv"),
            row.names=FALSE)
  write.csv(summary_df[, c("condition","total_weight")],
            file.path(out_subdir, "cellchat_interaction_weights.csv"),
            row.names=FALSE)

  # ── Merge + rankNet info-flow per pathway ───────────────────────────────
  cc_merge <- tryCatch(
    mergeCellChat(cc_list, add.names=names(cc_list)),
    error = function(e) {
      message("  mergeCellChat failed: ", conditionMessage(e))
      NULL
    })
  if (!is.null(cc_merge)) {
    saveRDS(cc_merge, file.path(out_subdir, "cellchat_merged.rds"))

    # rankNet returns a ggplot whose $data carries the per-pathway info-flow
    flow <- tryCatch(
      rankNet(cc_merge,
              mode    = "comparison",
              measure = "weight",
              comparison = seq_along(cc_list),
              stacked = FALSE, do.stat = TRUE, return.data = TRUE),
      error = function(e) {
        message("  rankNet failed: ", conditionMessage(e))
        NULL
      })
    if (!is.null(flow)) {
      flow_df <- if (is.list(flow) && !is.null(flow$signaling.contribution)) {
        flow$signaling.contribution
      } else if (inherits(flow, "ggplot")) {
        flow$data
      } else flow
      write.csv(flow_df,
                file.path(out_subdir, "cellchat_pathway_infoflow.csv"),
                row.names=FALSE)
    }
  }

  # ── Pairwise delta_interaction_counts (sender × receiver) ───────────────
  pair_df_all <- list()
  pairs <- list(
    c("Mutant_STM",  "Mutant_DMSO"),
    c("WT_STM",      "WT_DMSO"),
    c("Mutant_DMSO", "WT_DMSO"),
    c("Mutant_STM",  "WT_STM")
  )
  for (pr in pairs) {
    A <- pr[1]; B <- pr[2]
    if (!all(c(A,B) %in% names(cc_list))) next
    mA <- cc_list[[A]]@net$count
    mB <- cc_list[[B]]@net$count
    ct_common <- intersect(rownames(mA), rownames(mB))
    if (length(ct_common) < 2) next
    delta <- mA[ct_common, ct_common] - mB[ct_common, ct_common]
    d <- as.data.frame(as.table(delta))
    colnames(d) <- c("sender","receiver","delta_count")
    d$comparison <- paste0(A, "_vs_", B)
    pair_df_all[[paste0(A,"_vs_",B)]] <- d
  }
  if (length(pair_df_all) > 0) {
    pair_df <- do.call(rbind, pair_df_all)
    pair_df <- pair_df[order(pair_df$comparison, -abs(pair_df$delta_count)), ]
    write.csv(pair_df,
              file.path(out_subdir, "cellchat_delta_interaction_counts.csv"),
              row.names=FALSE)
  }

  # ── Differential L-R via subsetCommunication on merged ──────────────────
  if (!is.null(cc_merge)) {
    for (pr in pairs) {
      A <- pr[1]; B <- pr[2]
      if (!all(c(A,B) %in% names(cc_list))) next
      pair_name <- paste0(A, "_vs_", B)
      diff_lr <- tryCatch({
        netA <- subsetCommunication(cc_list[[A]])
        netB <- subsetCommunication(cc_list[[B]])
        key  <- c("source","target","interaction_name")
        m <- merge(netA, netB, by=key, all=TRUE, suffixes=c(paste0(".",A), paste0(".",B)))
        probA <- m[[paste0("prob.", A)]]; probA[is.na(probA)] <- 0
        probB <- m[[paste0("prob.", B)]]; probB[is.na(probB)] <- 0
        m$delta_prob <- probA - probB
        m[order(-abs(m$delta_prob)), c(key, paste0("prob.", A), paste0("prob.", B), "delta_prob")]
      }, error = function(e) {
        message("  diff L-R failed for ", pair_name, ": ", conditionMessage(e))
        NULL
      })
      if (!is.null(diff_lr)) {
        write.csv(diff_lr,
                  file.path(out_subdir, paste0("cellchat_diff_LR_", pair_name, ".csv")),
                  row.names=FALSE)
      }
    }
  }
}

# ── Run for both annotation levels ─────────────────────────────────────────
for (level_col in c("manual_level1", "manual_level2")) {
  if (!(level_col %in% colnames(meta))) {
    message("Skipping ", level_col, " — not present in metadata")
    next
  }
  run_level(level_col, file.path(args$out, level_col))
}

message("\nScript 16 (CellChat) complete. Outputs in: ", args$out)
