#!/usr/bin/env python3
"""
make_fig2_composition.py — Report Figure 2 (replacement).

The previous composition figure used (a) a stacked bar of the four conditions —
which makes every cell type look equally present — and (b) a linear y-axis
dominated by Neutrophils, hiding order-of-magnitude shifts in rare stem
populations. This redraws it so the effects are visible:

  Panel A: per-sample cell-type proportions on a log10 axis, coloured by the
           four conditions. Rare populations (HSC, STHSC) become legible.
  Panel B: Mutant-vs-WT log2 fold-change per cell type (pooled over treatment),
           sorted; cell types flagged credible by scCODA in the genotype
           contrast are filled, others are hollow.

Source: results/18_composition/final_annotation/ (canonical `cell_type`).

Usage:
  pixi run python report/make_fig2_composition.py
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D = "results/18_composition/final_annotation"
OUT = "report/figures/fig2_composition.png"
META = ["sample", "donor", "condition", "genotype", "treatment"]

COND_COLORS = {
    "Mutant_DMSO": "#1f77b4", "Mutant_STM": "#ff7f0e",
    "WT_DMSO": "#2ca02c", "WT_STM": "#d62728",
}


def genotype_credible_set():
    """Cell types credible for the genotype covariate in either treatment stratum."""
    cred = set()
    for f in glob.glob(os.path.join(D, "sccoda_genotype_in_*_credible_effects.csv")):
        df = pd.read_csv(f)
        g = df[df["covariate"].str.contains("genotype", case=False) & df["credible"]]
        cred |= set(g["celltype"])
    return cred


def main():
    c = pd.read_csv(os.path.join(D, "sccoda_counts_per_sample.csv"))
    ct = [col for col in c.columns if col not in META]
    prop = c[ct].div(c[ct].sum(1), axis=0)               # per-sample proportion
    cond = (c["genotype"].astype(str) + "_" + c["treatment"].astype(str)).values

    eps = 1e-4                                            # pseudocount for log axis
    # order cell types by overall abundance (most abundant first) for panel A
    order = prop.mean(0).sort_values(ascending=False).index.tolist()

    cred = genotype_credible_set()

    # Mutant-vs-WT log2FC, pooled over treatment
    mut = prop[c["genotype"] == "Mutant"].mean(0)
    wt = prop[c["genotype"] == "WT"].mean(0)
    lfc = np.log2((mut + eps) / (wt + eps)).reindex(order)

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(15, 7), gridspec_kw={"width_ratios": [2.2, 1]}
    )

    # ── Panel A: log-scale proportions, points per sample ───────────────
    x = np.arange(len(order))
    seen = set()
    for j, cn in enumerate(cond):
        col = COND_COLORS.get(cn, "gray")
        jitter = (j % 2 - 0.5) * 0.18 + np.random.uniform(-0.05, 0.05, len(order))
        y = (prop.iloc[j][order].values) + eps
        lbl = cn if cn not in seen else None
        seen.add(cn)
        axA.scatter(x + jitter, y, s=22, color=col, label=lbl,
                    edgecolors="none", alpha=0.85)
    axA.set_yscale("log")
    axA.set_xticks(x)
    axA.set_xticklabels(order, rotation=90, fontsize=7)
    axA.set_ylabel("Proportion of cells (log scale)", fontsize=9)
    axA.set_title("A. Cell-type proportions per sample", fontsize=11, loc="left")
    axA.legend(frameon=False, fontsize=8, markerscale=1.3, loc="upper right")
    axA.grid(axis="y", ls=":", alpha=0.4)

    # ── Panel B: Mutant vs WT log2FC lollipop ───────────────────────────
    lfc_sorted = lfc.sort_values()
    yb = np.arange(len(lfc_sorted))
    for yi, (name, val) in zip(yb, lfc_sorted.items()):
        is_cred = name in cred
        color = "#b2182b" if val > 0 else "#2166ac"
        axB.plot([0, val], [yi, yi], color=color, lw=1.2, zorder=1)
        axB.scatter(val, yi, s=46 if is_cred else 30,
                    facecolor=color if is_cred else "white",
                    edgecolor=color, linewidths=1.3, zorder=2)
    axB.axvline(0, color="black", lw=0.8)
    axB.set_yticks(yb)
    axB.set_yticklabels(lfc_sorted.index, fontsize=7)
    axB.set_xlabel("log2 fold-change  (Mutant / WT)", fontsize=9)
    axB.set_title("B. Genotype effect on composition", fontsize=11, loc="left")
    axB.grid(axis="x", ls=":", alpha=0.4)
    # legend for credibility
    from matplotlib.lines import Line2D
    axB.legend(handles=[
        Line2D([0], [0], marker="o", color="gray", markerfacecolor="gray",
               linestyle="none", markersize=7, label="scCODA-credible"),
        Line2D([0], [0], marker="o", color="gray", markerfacecolor="white",
               linestyle="none", markersize=6, label="not credible"),
    ], frameon=False, fontsize=7, loc="lower right")

    fig.suptitle("Cell-type composition (canonical `cell_type`, n=2 donors per group)",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"✓ wrote {OUT}")
    print("  genotype-credible:", sorted(cred))
    print("  top |log2FC| Mutant/WT:")
    print(lfc.reindex(lfc.abs().sort_values(ascending=False).index).head(8).round(2).to_string())


if __name__ == "__main__":
    np.random.seed(0)
    main()
