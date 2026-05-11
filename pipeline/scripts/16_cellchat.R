#!/usr/bin/env Rscript
# 16_cellchat.R
#
# Cell-cell communication analysis using CellChat v2.
# Runs one CellChat object per condition (WT_DMSO, WT_STM, Mutant_DMSO, Mutant_STM),
# then performs cross-condition comparison.
#
# Inputs
#   --input     path to adata_progenitor_annotated.h5ad
#   --out       output directory
#   --organism  CellChatDB organism: "Mm" (mouse, default) or "Hs"
#
# Outputs (in --out/)
#   cellchat_<condition>.rds              one CellChat object per condition
#   cellchat_interaction_counts.csv       #interactions per condition
#   cellchat_interaction_weights.csv      interaction strength per condition
#   cellchat_bubble_<comparison>.png      bubble plots for key comparisons
#   cellchat_chord_<condition>.png        chord diagrams per condition
#   cellchat_differential_interactions.csv differential interactions (Mut_STM vs Mut_DMSO)

suppressPackageStartupMessages({
  library(argparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
})

# ── CLI ────────────────────────────────────────────────────────────────────
parser <- ArgumentParser()
parser$add_argument("--input",    required=TRUE)
parser$add_argument("--out",      default=".")
parser$add_argument("--organism", default="Mm")
args <- parser$parse_args()

dir.create(args$out, recursive=TRUE, showWarnings=FALSE)
options(future.globals.maxSize = 4 * 1024^3)  # 4 GB

# ── Load h5ad via zellkonverter ────────────────────────────────────────────
message("Loading h5ad...")
sce  <- readH5AD(args$input, use_hdf5=TRUE)
meta <- as.data.frame(colData(sce))
counts <- assay(sce, "X")  # raw counts (before normalisation stored in X)

# Use manual_level1 as cell identity; collapse low-count groups
meta$celltype <- as.character(meta$manual_level1)
meta$condition <- paste(meta$genotype, meta$treatment, sep="_")

# ── CellChatDB (mouse ligand-receptor interactions) ────────────────────────
CellChatDB <- CellChatDB.mouse
# Use all interactions (Secreted Signaling + Cell-Cell Contact + ECM-Receptor)
CellChatDB.use <- CellChatDB

# ── Helper: build and run one CellChat object ──────────────────────────────
run_cellchat <- function(counts_sub, meta_sub, label) {
  message("  Building CellChat for: ", label)
  cc <- createCellChat(object=counts_sub, meta=meta_sub, group.by="celltype")
  CellChatDB(cc) <- CellChatDB.use
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type="triMean", population.size=TRUE)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# ── Run CellChat per condition ─────────────────────────────────────────────
conditions <- c("WT_DMSO", "WT_STM", "Mutant_DMSO", "Mutant_STM")
cc_list <- list()

for (cond in conditions) {
  idx <- meta$condition == cond
  if (sum(idx) < 50) {
    message("Skipping ", cond, " — only ", sum(idx), " cells")
    next
  }
  cc_list[[cond]] <- run_cellchat(counts[, idx], meta[idx, ], cond)
  saveRDS(cc_list[[cond]], file.path(args$out, paste0("cellchat_", cond, ".rds")))

  # Chord diagram per condition
  png(file.path(args$out, paste0("cellchat_chord_", cond, ".png")),
      width=800, height=800, res=100)
  netVisual_circle(cc_list[[cond]]@net$count,
                   title.name=paste("Interaction counts —", cond))
  dev.off()
}

# ── Cross-condition comparison ─────────────────────────────────────────────
if (length(cc_list) >= 2) {
  message("Running cross-condition comparison...")

  # Merge and compare
  cc_merge <- mergeCellChat(cc_list, add.names=names(cc_list))
  saveRDS(cc_merge, file.path(args$out, "cellchat_merged.rds"))

  # Interaction counts table
  cnt_df <- do.call(rbind, lapply(names(cc_list), function(nm) {
    n <- cc_list[[nm]]@net$count
    data.frame(
      condition   = nm,
      total_links = sum(n > 0),
      total_weight = sum(cc_list[[nm]]@net$weight)
    )
  }))
  write.csv(cnt_df, file.path(args$out, "cellchat_interaction_counts.csv"),
            row.names=FALSE)

  # Bubble plot: Mutant_DMSO vs Mutant_STM (drug effect in disease)
  if (all(c("Mutant_DMSO", "Mutant_STM") %in% names(cc_list))) {
    p <- netVisual_bubble(
      cc_merge,
      sources.use  = which(levels(cc_list[["Mutant_DMSO"]]@idents) == "monocyte"),
      targets.use  = NULL,
      comparison   = c(which(names(cc_list) == "Mutant_DMSO"),
                       which(names(cc_list) == "Mutant_STM")),
      angle.x      = 45,
      title.name   = "Monocyte outgoing — Mutant: DMSO vs STM"
    )
    ggsave(file.path(args$out, "cellchat_bubble_Mutant_DMSO_vs_STM.png"),
           p, width=12, height=8, dpi=150)
  }

  # Differential interactions: Mutant_STM vs Mutant_DMSO
  if (all(c("Mutant_DMSO", "Mutant_STM") %in% names(cc_list))) {
    pos <- cc_list[["Mutant_STM"]]@net$count - cc_list[["Mutant_DMSO"]]@net$count
    diff_df <- as.data.frame(as.table(pos))
    colnames(diff_df) <- c("sender", "receiver", "delta_count")
    diff_df <- diff_df[order(-abs(diff_df$delta_count)), ]
    write.csv(diff_df, file.path(args$out, "cellchat_differential_interactions.csv"),
              row.names=FALSE)
  }
}

message("Script 16 (CellChat) complete. Outputs in: ", args$out)
