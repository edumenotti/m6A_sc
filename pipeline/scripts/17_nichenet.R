#!/usr/bin/env Rscript
# 17_nichenet.R
#
# NicheNet ligand activity analysis.
# Sender cells:   monocyte (manual_level1 == "monocyte")
# Receiver cells: progenitor + myeloid_progenitor (manual_level1)
#
# For each pairwise comparison, identifies top-ranked sender ligands predicted
# to regulate receiver DEGs, and scores ligand-target gene relationships.
#
# Inputs
#   --input       path to adata_progenitor_annotated.h5ad
#   --out         output directory
#   --organism    "mouse" (default) or "human"
#
# Outputs (in --out/)
#   nichenet_ligand_activity_<comparison>.csv   ligand activity scores
#   nichenet_ligand_target_heatmap_<comparison>.png
#   nichenet_receptor_expression_<comparison>.png

suppressPackageStartupMessages({
  library(argparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(nichenetr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

parser <- ArgumentParser()
parser$add_argument("--input",    required=TRUE)
parser$add_argument("--out",      default=".")
parser$add_argument("--organism", default="mouse")
args <- parser$parse_args()

dir.create(args$out, recursive=TRUE, showWarnings=FALSE)

# ── Load NicheNet prior networks ───────────────────────────────────────────
message("Loading NicheNet networks (mouse)...")
organism <- args$organism
if (organism == "mouse") {
  lr_network    <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_net  <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
} else {
  lr_network    <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_net  <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
}
lr_network <- lr_network %>% distinct(from, to)

# ── Load h5ad ──────────────────────────────────────────────────────────────
message("Loading h5ad...")
sce  <- readH5AD(args$input, use_hdf5=TRUE)
meta <- as.data.frame(colData(sce))
# normalised log1p counts expected in "X"; if raw use logcounts
expr <- as.matrix(assay(sce, "X"))
rownames(expr) <- rownames(sce)
colnames(expr) <- colnames(sce)
meta$condition <- paste(meta$genotype, meta$treatment, sep="_")

# ── Helper: run NicheNet for one comparison ────────────────────────────────
run_nichenet <- function(expr, meta, cond_sender, cond_receiver,
                         label, top_n_ligands=20) {
  message("  NicheNet: ", label)

  sender_idx   <- meta$manual_level1 == "monocyte" & meta$condition == cond_sender
  receiver_idx <- meta$manual_level1 %in% c("progenitor", "myeloid_progenitor") &
                  meta$condition == cond_receiver
  ctrl_idx     <- meta$manual_level1 %in% c("progenitor", "myeloid_progenitor") &
                  meta$condition != cond_receiver

  if (sum(sender_idx) < 10 || sum(receiver_idx) < 10) {
    message("    Skipping — insufficient cells (sender=",
            sum(sender_idx), ", receiver=", sum(receiver_idx), ")")
    return(NULL)
  }

  # Expressed genes per compartment
  sender_expressed   <- rownames(expr)[rowMeans(expr[, sender_idx, drop=FALSE] > 0) > 0.10]
  receiver_expressed <- rownames(expr)[rowMeans(expr[, receiver_idx, drop=FALSE] > 0) > 0.10]

  # DE genes in receiver (condition vs rest) using simple fold-change threshold
  rec_mean <- rowMeans(expr[, receiver_idx, drop=FALSE])
  ctrl_mean <- rowMeans(expr[, ctrl_idx, drop=FALSE]) + 0.001
  lfc <- log2((rec_mean + 0.001) / ctrl_mean)
  geneset_oi <- names(lfc[lfc > 0.5 & rec_mean > 0.05])
  geneset_oi <- intersect(geneset_oi, rownames(ligand_target))

  if (length(geneset_oi) < 10) {
    message("    Skipping — fewer than 10 DE genes in receiver")
    return(NULL)
  }

  # Ligands expressed by sender
  ligands_all  <- unique(lr_network$from)
  receptors_all <- unique(lr_network$to)
  expressed_ligands   <- intersect(sender_expressed, ligands_all)
  expressed_receptors <- intersect(receiver_expressed, receptors_all)

  potential_ligands <- lr_network %>%
    filter(from %in% expressed_ligands, to %in% expressed_receptors) %>%
    pull(from) %>% unique()

  if (length(potential_ligands) < 3) {
    message("    Skipping — fewer than 3 potential ligands")
    return(NULL)
  }

  # Ligand activity
  bg_genes <- rownames(ligand_target)
  activity <- predict_ligand_activities(
    geneset_oi         = geneset_oi,
    background_expressed_genes = intersect(receiver_expressed, bg_genes),
    ligand_target_matrix = ligand_target,
    potential_ligands  = potential_ligands
  )
  activity <- activity %>% arrange(desc(pearson))
  write.csv(activity,
            file.path(args$out, paste0("nichenet_ligand_activity_", label, ".csv")),
            row.names=FALSE)

  # Top ligand–target heatmap
  top_ligands <- head(activity$test_ligand, top_n_ligands)
  active_lt   <- ligand_target[top_ligands, geneset_oi, drop=FALSE]
  lt_df <- as.data.frame(active_lt) %>%
    tibble::rownames_to_column("ligand") %>%
    pivot_longer(-ligand, names_to="target", values_to="score") %>%
    filter(score > 0.001)

  p <- ggplot(lt_df, aes(target, ligand, fill=score)) +
    geom_tile() +
    scale_fill_gradient(low="white", high="#d62728") +
    theme_classic(base_size=8) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +
    labs(title=paste("Ligand-target matrix —", label),
         x="Target gene (receiver DE)", y="Ligand (sender)")
  ggsave(file.path(args$out,
                   paste0("nichenet_ligand_target_heatmap_", label, ".png")),
         p, width=12, height=6, dpi=150)

  activity
}

# ── Run comparisons ────────────────────────────────────────────────────────
# Primary: monocytes in Mutant_STM → progenitors in Mutant_STM
#          (what are monocytes signaling to MDS progenitors during treatment?)
run_nichenet(expr, meta,
             cond_sender   = "Mutant_STM",
             cond_receiver = "Mutant_STM",
             label         = "Mutant_STM_mono_to_prog")

# Reference: WT_STM monocytes → WT_STM progenitors
run_nichenet(expr, meta,
             cond_sender   = "WT_STM",
             cond_receiver = "WT_STM",
             label         = "WT_STM_mono_to_prog")

# Baseline disease: Mutant_DMSO monocytes → Mutant_DMSO progenitors
run_nichenet(expr, meta,
             cond_sender   = "Mutant_DMSO",
             cond_receiver = "Mutant_DMSO",
             label         = "Mutant_DMSO_mono_to_prog")

message("Script 17 (NicheNet) complete. Outputs in: ", args$out)
