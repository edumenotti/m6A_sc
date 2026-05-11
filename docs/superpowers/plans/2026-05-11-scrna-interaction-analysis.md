# scRNA-seq Interaction Analysis (Scripts 15–17) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add three downstream analysis scripts (macrophage states, CellChat, NicheNet) that address the core biological question — what is the role of macrophages in METTL3-inhibitor-driven clearance of MDS cells — and wire them into the Nextflow pipeline as an optional `run_interaction_analysis` block.

**Architecture:** All three analyses take `adata_progenitor_annotated.h5ad` (with `genotype` ∈ {WT, Mutant} and `treatment` ∈ {DMSO, STM}) as input and run in parallel after `APPLY_PROGENITOR_ANNOTATION`. Scripts 15 is Python (uses existing scanpy/decoupler/pydeseq2); scripts 16 and 17 are R (use zellkonverter already in pixi, plus CellChat and NicheNet installed via GitHub in Apptainer). The block is gated by `params.run_interaction_analysis` (default `false`).

**Tech Stack:** Python 3.11 + scanpy + decoupler + pydeseq2; R 4.5 + zellkonverter + CellChat v2 + NicheNet (nichenetr); Nextflow DSL2; Apptainer (ubuntu:22.04 + pixi); pixi (conda-forge + bioconda).

---

## Biological Context

The four experimental conditions in the dataset are:

| Condition | Meaning |
|---|---|
| `WT_DMSO` | Wild-type cells, no treatment (baseline) |
| `WT_STM` | Wild-type cells + METTL3 inhibitor |
| `Mutant_DMSO` | SRSF2/IDH2-mutant MDS cells, no treatment |
| `Mutant_STM` | SRSF2/IDH2-mutant MDS cells + METTL3 inhibitor |

The central comparison is **Mutant_STM vs Mutant_DMSO** (drug effect in disease), with **WT_STM vs WT_DMSO** as a biological reference. The column `genotype` was derived from `replicate` (1=WT/CD45.1, 2=Mutant/CD45.2) and is present in the h5ad.

Cell type labels live in `manual_level1` and `manual_level2`. The monocyte cluster (`manual_level1 == "monocyte"`, `manual_level2 == "Classical_monocyte"`) is the proxy for macrophage/myeloid innate immune cells in this dataset (bone marrow macrophages are not well represented in CD45+ sort; Classical_monocyte = bone marrow-resident monocyte / macrophage precursor). Progenitor compartment = `manual_level1` ∈ {`progenitor`, `myeloid_progenitor`, `cycling_myeloid`}.

---

## File Map

**New scripts:**
- `pipeline/scripts/15_macrophage_states.py` — pseudobulk DEG + M1/M2/OAS scoring in monocyte cluster
- `pipeline/scripts/16_cellchat.R` — CellChat v2 cell-cell communication across 4 conditions
- `pipeline/scripts/17_nichenet.R` — NicheNet ligand activity (monocyte → progenitor signaling)

**New Nextflow modules:**
- `pipeline/modules/macrophage_states.nf`
- `pipeline/modules/cellchat.nf`
- `pipeline/modules/nichenet.nf`

**Modified files:**
- `pipeline/main.nf` — add `run_interaction_analysis` block
- `pipeline/params.yaml` — add `run_interaction_analysis: false`, `cellchat_organism: "Mm"`, `nichenet_organism: "mouse"`
- `pipeline/nextflow.config` — add CELLCHAT + NICHENET to general partition in slurm profile
- `pixi.toml` — add `bioconductor-muscat`
- `charles-scrna.def` — add CellChat v2 and NicheNet GitHub installs + smoke tests

---

## Task 1: Verify R dependencies in pixi.toml

Scripts 16 and 17 use R packages already present in `pixi.toml` (`zellkonverter`, `r-seurat`, `bioconductor-singlecellexperiment`). CellChat and NicheNet are installed via GitHub in the Apptainer def (Task 2). No new pixi dependencies are needed — this task just confirms the environment is consistent.

**Files:**
- Read: `pixi.toml`

- [ ] **Step 1: Confirm required R packages are present**

```bash
grep -E "zellkonverter|r-seurat|singlecellexperiment|r-base" /home/edu-pc/Yale/Charles/pixi.toml
```

Expected output includes all four. If any are missing, add them under `[dependencies]` in `pixi.toml` and re-run `pixi install`.

- [ ] **Step 2: Smoke test R can load the data bridge**

```bash
pixi run Rscript -e "library(zellkonverter); library(SingleCellExperiment); cat('R bridge ok\n')"
```

Expected output: `R bridge ok`

- [ ] **Step 3: Commit if pixi.toml was modified**

```bash
# Only if pixi.toml was changed in Step 1:
git add pixi.toml pixi.lock
git commit -m "deps: ensure R bridge packages present for interaction analysis"
```

---

## Task 2: Update Apptainer def — install CellChat v2 and NicheNet

CellChat and NicheNet are not on conda-forge/bioconda; they must be installed from GitHub inside the container.

**Files:**
- Modify: `charles-scrna.def`

- [ ] **Step 1: Add GitHub R installs to the `%post` section**

Open `charles-scrna.def`. Locate the R package installation block (after `pixi install --locked`). Add the following block immediately after it:

```bash
    echo "=== Installing CellChat v2 and NicheNet from GitHub ===" && \
    Rscript -e "
      options(repos = c(CRAN='https://cloud.r-project.org'))
      if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes')
      remotes::install_github('jinworks/CellChat', upgrade='never', quiet=TRUE)
      remotes::install_github('saeyslab/nichenetr', upgrade='never', quiet=TRUE)
      cat('CellChat + NicheNet installed\n')
    " && \
```

- [ ] **Step 2: Add smoke tests to the `%test` section**

Locate the `%test` section. After the existing `Rscript` smoke test line, add:

```bash
Rscript -e "library(CellChat); library(nichenetr); cat('CellChat + NicheNet ok\n')"
```

- [ ] **Step 3: Commit the def file (container rebuild happens on HPC)**

```bash
git add charles-scrna.def
git commit -m "container: add CellChat v2 and NicheNet to Apptainer def"
```

> **Note:** The `.sif` is in `.gitignore`. Rebuild on HPC login node with:
> `apptainer build --fakeroot charles-scrna.sif charles-scrna.def`
> This requires internet access — run on a login node, not a compute node.

---

## Task 3: Script 15 — Monocyte state scoring (Python)

Answers: *Are monocytes collectively shifting toward M1/OAS state with STM treatment? Is this shift blunted in the Mutant genotype vs WT?*

**Approach:** Gene set scoring with `decoupler` ULM on the monocyte cluster (~1096 cells). No pseudobulk DEG — with only 2 biological replicates (D1, D2) the statistical power is insufficient for formal testing. The scoring is single-cell-level and produces interpretable distributions per condition.

**Files:**
- Create: `pipeline/scripts/15_macrophage_states.py`

- [ ] **Step 1: Create the script**

```python
#!/usr/bin/env python3
"""
15_macrophage_states.py

Gene-set scoring for the monocyte cluster (manual_level1 == "monocyte"),
stratified across 4 conditions: WT_DMSO, WT_STM, Mutant_DMSO, Mutant_STM.

Uses decoupler ULM to score M1 (pro-inflammatory), M2 (anti-inflammatory),
and OAS/interferon programs per cell, then plots distributions per condition.

Note: pseudobulk DEG is intentionally omitted — with n=2 donors (D1, D2)
the statistical power is insufficient for reliable formal testing.

Inputs
------
--input  : path to adata_progenitor_annotated.h5ad (must have genotype column)
--out    : output directory

Outputs
-------
monocyte_state_scores.csv   : per-cell scores (M1, M2, OAS) + condition label
monocyte_state_scores.png   : boxplots of each score across 4 conditions
monocyte_score_summary.csv  : median score per condition (for reporting)
"""
import argparse
import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.stats as stats

warnings.filterwarnings("ignore", category=FutureWarning)
sc.settings.verbosity = 1

# ── Gene sets (mouse symbols) ──────────────────────────────────────────────
# M1: classical pro-inflammatory / anti-tumoral activation
M1_GENES = [
    "Nos2", "Tnf", "Il1b", "Il6", "Il12b", "Cxcl9", "Cxcl10",
    "Cd80", "Cd86", "Stat1", "Irf5", "H2-Aa", "H2-Ab1",
]
# M2: anti-inflammatory / immunosuppressive / pro-tumoral
M2_GENES = [
    "Arg1", "Mrc1", "Cd163", "Il10", "Tgfb1",
    "Mgl2", "Ccl22", "Stat6", "Irf4", "Pparg",
]
# OAS/interferon: dsRNA sensing downstream of METTL3 inhibition
OAS_GENES = [
    "Oas1a", "Oas2", "Oas3", "Oasl1", "Oasl2",
    "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Isg15",
]

CONDITIONS = ["WT_DMSO", "WT_STM", "Mutant_DMSO", "Mutant_STM"]
COLORS     = ["#4c72b0", "#55a868", "#c44e52", "#dd8452"]


def score_gene_sets(adata: sc.AnnData, gene_sets: dict) -> pd.DataFrame:
    """Score gene sets per cell using decoupler ULM. Returns DataFrame of scores."""
    net = pd.concat([
        pd.DataFrame({"source": name, "target": genes, "weight": 1.0})
        for name, genes in gene_sets.items()
    ])
    net = net[net["target"].isin(adata.var_names)]
    missing = {
        name: [g for g in genes if g not in adata.var_names]
        for name, genes in gene_sets.items()
    }
    for name, absent in missing.items():
        if absent:
            print(f"  [{name}] genes not in dataset: {absent}")

    dc.run_ulm(adata, net=net, source="source", target="target",
               weight="weight", use_raw=False, verbose=False)
    score_cols = list(gene_sets.keys())
    scores = adata.obsm["ulm_estimate"][score_cols].copy()
    scores.index = adata.obs.index
    return scores


def plot_scores(scores: pd.DataFrame, out_path: str) -> None:
    """Boxplot of M1/M2/OAS scores per condition."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    for ax, sig_name in zip(axes, ["M1", "M2", "OAS"]):
        data = [scores.loc[scores["condition"] == c, sig_name].dropna().values
                for c in CONDITIONS]
        bp = ax.boxplot(data, patch_artist=True, notch=False, widths=0.55,
                        medianprops=dict(color="black", linewidth=2),
                        flierprops=dict(marker=".", markersize=2, alpha=0.3))
        for patch, color in zip(bp["boxes"], COLORS):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        ax.set_xticks(range(1, len(CONDITIONS) + 1))
        ax.set_xticklabels(CONDITIONS, rotation=25, ha="right", fontsize=8)
        ax.set_title(f"{sig_name} program score")
        ax.set_ylabel("decoupler ULM score")

        # Annotate with median per condition
        for i, d in enumerate(data, start=1):
            if len(d):
                ax.text(i, np.median(d), f"{np.median(d):.2f}",
                        ha="center", va="bottom", fontsize=7, color="black")

    fig.suptitle("Monocyte transcriptional state — 4 conditions", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out",   default=".")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ── Subset to monocyte cluster ──────────────────────────────────────
    mac = adata[adata.obs["manual_level1"] == "monocyte"].copy()
    mac.obs["condition"] = (
        mac.obs["genotype"].astype(str) + "_" + mac.obs["treatment"].astype(str)
    )
    print(f"Monocyte cells: {mac.shape[0]}")
    print(mac.obs["condition"].value_counts())

    # ── Score gene sets ─────────────────────────────────────────────────
    scores = score_gene_sets(mac, {"M1": M1_GENES, "M2": M2_GENES, "OAS": OAS_GENES})
    scores["condition"] = mac.obs["condition"].values
    scores.to_csv(os.path.join(args.out, "monocyte_state_scores.csv"))

    # ── Summary table: median per condition ────────────────────────────
    summary = scores.groupby("condition")[["M1", "M2", "OAS"]].median()
    summary.index.name = "condition"
    summary.to_csv(os.path.join(args.out, "monocyte_score_summary.csv"))
    print("\nMedian scores per condition:")
    print(summary.round(3).to_string())

    # ── Plots ───────────────────────────────────────────────────────────
    plot_scores(scores, os.path.join(args.out, "monocyte_state_scores.png"))

    print("Script 15 complete.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke test locally**

```bash
cd /home/edu-pc/Yale/Charles
pixi run python pipeline/scripts/15_macrophage_states.py \
    --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
    --out results/15_macrophage_states
```

Expected: creates `results/15_macrophage_states/` with 3 files. Verify:

```bash
ls results/15_macrophage_states/
# monocyte_state_scores.csv  monocyte_state_scores.png  monocyte_score_summary.csv
cat results/15_macrophage_states/monocyte_score_summary.csv
```

Expected summary shows M1 and OAS scores higher in `WT_STM` and `Mutant_STM` vs their DMSO counterparts, with the shift potentially smaller in the Mutant group (consistent with the blunted OAS response described in the meeting).

- [ ] **Step 3: Update Nextflow module output list to match new filenames**

In `pipeline/modules/macrophage_states.nf` (created in Task 6), the output block references `macrophage_state_scores.csv`. Update it to:

```groovy
    output:
    path "monocyte_state_scores.csv"
    path "monocyte_score_summary.csv"
    path "*.png"
```

- [ ] **Step 4: Commit**

```bash
git add pipeline/scripts/15_macrophage_states.py
git commit -m "feat(15): monocyte state scoring with decoupler ULM (no pseudobulk)"
```

---

## Task 4: Script 16 — CellChat (R)

Answers: *Which cell types are talking to which, and does the STM treatment change this? Are monocytes sending or receiving more signals in MDS vs WT?*

**Files:**
- Create: `pipeline/scripts/16_cellchat.R`

- [ ] **Step 1: Create the script**

```r
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
```

- [ ] **Step 2: Smoke test locally (requires R + CellChat installed)**

```bash
pixi run Rscript pipeline/scripts/16_cellchat.R \
    --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
    --out results/16_cellchat
```

Expected: runs without error and creates `results/16_cellchat/` with at minimum `cellchat_WT_DMSO.rds` and `cellchat_interaction_counts.csv`. If CellChat is not yet installed locally (only in Apptainer), this step runs on the HPC after the container is rebuilt (Task 2). Verify by running in the container:

```bash
apptainer exec charles-scrna.sif Rscript pipeline/scripts/16_cellchat.R \
    --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
    --out results/16_cellchat
```

- [ ] **Step 3: Commit**

```bash
git add pipeline/scripts/16_cellchat.R
git commit -m "feat(16): CellChat cell-cell communication analysis"
```

---

## Task 5: Script 17 — NicheNet (R)

Answers: *Which specific ligands from monocytes are predicted to regulate gene expression in MDS progenitors? How does this change with STM treatment?*

**Files:**
- Create: `pipeline/scripts/17_nichenet.R`

- [ ] **Step 1: Create the script**

```r
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
```

- [ ] **Step 2: Smoke test (requires NicheNet + internet for network files)**

```bash
pixi run Rscript pipeline/scripts/17_nichenet.R \
    --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
    --out results/17_nichenet
```

If running locally without NicheNet, test on HPC after container rebuild:

```bash
apptainer exec charles-scrna.sif Rscript pipeline/scripts/17_nichenet.R \
    --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
    --out results/17_nichenet
```

Expected: `results/17_nichenet/` contains at minimum `nichenet_ligand_activity_Mutant_STM_mono_to_prog.csv` and the corresponding heatmap PNG.

> **HPC note:** NicheNet downloads network files (~150 MB) from Zenodo on first run. On Grace, this requires internet access on the login node. If compute nodes lack internet, pre-download the 3 RDS files and add a `--networks_dir` argument pointing to a local cache.

- [ ] **Step 3: Commit**

```bash
git add pipeline/scripts/17_nichenet.R
git commit -m "feat(17): NicheNet ligand-receptor activity analysis"
```

---

## Task 6: Nextflow modules for scripts 15–17

**Files:**
- Create: `pipeline/modules/macrophage_states.nf`
- Create: `pipeline/modules/cellchat.nf`
- Create: `pipeline/modules/nichenet.nf`

- [ ] **Step 1: Create `macrophage_states.nf`**

```groovy
process MACROPHAGE_STATES {
    tag "macrophage_states"
    publishDir "${params.outdir}/15_macrophage_states", mode: 'copy'
    memory '16 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "macrophage_state_scores.csv"
    path "macrophage_pseudobulk_deg.csv"
    path "*.png"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/15_macrophage_states.py \
        --input ${h5ad} \
        --out .
    """
}
```

- [ ] **Step 2: Create `cellchat.nf`**

```groovy
process CELLCHAT {
    tag "cellchat"
    publishDir "${params.outdir}/16_cellchat", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "cellchat_*.rds"
    path "cellchat_interaction_counts.csv"
    path "cellchat_differential_interactions.csv", optional: true
    path "*.png"

    script:
    """
    pixi run -m ${params.pixi_manifest} Rscript ${projectDir}/scripts/16_cellchat.R \
        --input ${h5ad} \
        --out . \
        --organism ${params.cellchat_organism}
    """
}
```

- [ ] **Step 3: Create `nichenet.nf`**

```groovy
process NICHENET {
    tag "nichenet"
    publishDir "${params.outdir}/17_nichenet", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "nichenet_ligand_activity_*.csv"
    path "*.png"

    script:
    """
    pixi run -m ${params.pixi_manifest} Rscript ${projectDir}/scripts/17_nichenet.R \
        --input ${h5ad} \
        --out . \
        --organism ${params.nichenet_organism}
    """
}
```

- [ ] **Step 4: Commit**

```bash
git add pipeline/modules/macrophage_states.nf \
        pipeline/modules/cellchat.nf \
        pipeline/modules/nichenet.nf
git commit -m "feat: Nextflow modules for scripts 15-17 (interaction analysis)"
```

---

## Task 7: Wire into main.nf and update params.yaml

**Files:**
- Modify: `pipeline/main.nf`
- Modify: `pipeline/params.yaml`
- Modify: `pipeline/nextflow.config`

- [ ] **Step 1: Add params to `params.yaml`**

Append to the end of `pipeline/params.yaml`:

```yaml
# ── Interaction analysis (scripts 15–17) ──
run_interaction_analysis: false
cellchat_organism: "Mm"
nichenet_organism: "mouse"
```

- [ ] **Step 2: Add includes and conditional block to `main.nf`**

Add to the `include` block at the top of `main.nf` (after the existing includes):

```groovy
include { MACROPHAGE_STATES } from './modules/macrophage_states'
include { CELLCHAT }          from './modules/cellchat'
include { NICHENET }          from './modules/nichenet'
```

Add to the end of the `workflow {}` block, after the progenitor recluster block:

```groovy
    /*
     * Optional interaction analysis block (run_interaction_analysis=true):
     *   Requires APPLY_PROGENITOR_ANNOTATION to have produced the annotated h5ad
     *   with genotype column (replicate 1=WT, 2=Mutant).
     *   MACROPHAGE_STATES, CELLCHAT, and NICHENET run in parallel.
     */
    if (params.run_interaction_analysis) {
        map_file = file(params.progenitor_annotation_map)
        if (!map_file.exists()) {
            error "run_interaction_analysis=true requires progenitor_annotation_map to exist. Run the progenitor sub-workflow first."
        }
        prog_annotated_ch = Channel.fromPath(
            "${params.outdir}/14_progenitor_annotated/adata_progenitor_annotated.h5ad",
            checkIfExists: true
        )
        MACROPHAGE_STATES(prog_annotated_ch)
        CELLCHAT(prog_annotated_ch)
        NICHENET(prog_annotated_ch)
    }
```

- [ ] **Step 3: Add CELLCHAT and NICHENET to slurm general partition in `nextflow.config`**

Open `pipeline/nextflow.config`. Locate the slurm profile `withLabel: 'general'` block and add:

```groovy
        withName: 'MACROPHAGE_STATES|CELLCHAT|NICHENET' {
            clusterOptions = '--partition=general --time=04:00:00'
        }
```

- [ ] **Step 4: Validate DAG with stub run**

```bash
cd /home/edu-pc/Yale/Charles/pipeline
nextflow run main.nf -profile local -params-file params.yaml \
    --run_interaction_analysis true -stub 2>&1 | grep -E "MACRO|CELL|NICHE|ERROR|error"
```

Expected: sees `MACROPHAGE_STATES`, `CELLCHAT`, `NICHENET` in the stub output with no errors.

- [ ] **Step 5: Commit**

```bash
git add pipeline/main.nf pipeline/params.yaml pipeline/nextflow.config
git commit -m "feat: wire interaction analysis block into main.nf (scripts 15-17)"
```

---

## Task 8: Final integration commit and GitHub push readiness

- [ ] **Step 1: Verify git status is clean**

```bash
git status
git log --oneline -8
```

Expected: no untracked `.rds`, `.h5ad`, or `*.png` from results/ (these should be in `.gitignore`). Only pipeline scripts, modules, config, and def file are tracked.

- [ ] **Step 2: Verify `.gitignore` covers results and container**

```bash
grep -E "^results/|^.*\.sif$|^.*\.rds$" .gitignore
```

Expected output includes lines for `results/`, `*.sif`, and optionally `*.rds`. If `*.rds` is missing (CellChat objects can be large), add it:

```bash
echo "*.rds" >> .gitignore
git add .gitignore
git commit -m "chore: ignore *.rds files (CellChat objects)"
```

- [ ] **Step 3: Push to GitHub**

```bash
git push origin main
```

---

## Self-Review

### Spec coverage

| Requirement | Covered by |
|---|---|
| CellChat interaction analysis | Tasks 4, 6, 7 |
| Macrophage state analysis | Tasks 3, 6, 7 |
| NicheNet mechanistic follow-up | Tasks 5, 6, 7 |
| Pipeline integration (Nextflow) | Tasks 6, 7 |
| Apptainer compatibility | Task 2 |
| GitHub record | Task 8 |
| `run_interaction_analysis` flag | Task 7 |
| 4 conditions (WT/Mutant × DMSO/STM) | Scripts 15, 16, 17 all use `genotype` + `treatment` |

### Placeholder scan

No TBD, TODO, or "implement later" in any code block. All file paths are absolute or use `${projectDir}` as in the rest of the pipeline.

### Type consistency

- `manual_level1` used consistently across all 3 scripts to filter monocytes and progenitors
- `genotype` column (added in pre-analysis step) used in scripts 15 and 16; script 17 derives `condition` = `genotype + "_" + treatment` identically to script 16
- `params.cellchat_organism` passed as `--organism` to script 16; `params.nichenet_organism` to script 17
- All Nextflow modules use `${params.pixi_manifest}` and `${projectDir}/scripts/` consistent with existing modules
