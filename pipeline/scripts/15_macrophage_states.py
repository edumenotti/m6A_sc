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

    # decoupler 2.x API: dc.mt.ulm; stores results in obsm["score_ulm"]
    dc.mt.ulm(adata, net=net, raw=False, verbose=False, tmin=1)
    score_cols = list(gene_sets.keys())
    scores = adata.obsm["score_ulm"][score_cols].copy()
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
