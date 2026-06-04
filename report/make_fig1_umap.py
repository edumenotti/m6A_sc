#!/usr/bin/env python3
"""
make_fig1_umap.py — Report Figure 1.

Regenerate the UMAP landscape coloured by the CANONICAL `cell_type` annotation
(32 types; HemaScribe broad/fine + manual rescue, script 21z). The earlier
figure was coloured by manual_level1/level2 and so did not match the 32-type
count cited in the report.

Reads the embedding (obsm['X_umap']) straight from the final object; no
recomputation. Right-margin legend + on-data cluster labels for readability.

Usage:
  pixi run python report/make_fig1_umap.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import scanpy as sc

H5AD = "results/14_progenitor_annotated/adata_hemascribe.h5ad"
OUT = "report/figures/fig1_umap_celltypes.png"
KEY = "cell_type"

# Biologically ordered so the legend reads stem -> progenitor -> mature.
ORDER = [
    "HSC", "STHSC", "MPP2", "MPP4", "FcG_neg_MPP3", "FcG_pos_MPP3",
    "GMP", "GP", "mGMP", "cMoP", "Monocyte", "cDC", "pDC",
    "Immature_neutrophil", "Neutrophil", "Basophil_prog", "Basophil",
    "MkP", "Megakaryocyte", "MEP_Erythroid", "EryP", "RBC",
    "CLP", "Pro_Pre_B", "Immature_B_cell", "B_cell", "Plasma_cell",
    "T_cell", "NK_cell",
    "Progenitor_ambiguous", "Ambiguous", "NotHem",
]


def build_palette(categories):
    """32 visually distinct colours by concatenating tab20 + tab20b."""
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    return {c: base[i % len(base)] for i, c in enumerate(categories)}


def main():
    adata = sc.read_h5ad(H5AD)
    if "X_umap" not in adata.obsm:
        raise SystemExit("No X_umap embedding in the object.")

    # apply the biological order (keep any unexpected labels at the end)
    cats = [c for c in ORDER if c in set(adata.obs[KEY].astype(str))]
    cats += [c for c in adata.obs[KEY].astype(str).unique() if c not in cats]
    adata.obs[KEY] = adata.obs[KEY].astype(str).astype("category").cat.reorder_categories(cats)
    palette = build_palette(cats)

    xy = adata.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(11, 8))
    for c in cats:
        m = (adata.obs[KEY] == c).values
        ax.scatter(xy[m, 0], xy[m, 1], s=2.5, linewidths=0,
                   color=palette[c], label=c, rasterized=True)
        # on-data label at the cluster median, with a white halo for contrast
        ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), c,
                fontsize=5.5, ha="center", va="center", color="black",
                fontweight="bold",
                path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])

    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("UMAP1", fontsize=9); ax.set_ylabel("UMAP2", fontsize=9)
    ax.set_title(f"Cell-type landscape (canonical `cell_type`, {len(cats)} types)",
                 fontsize=11)
    leg = ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5),
                    frameon=False, fontsize=6.5, markerscale=3,
                    handletextpad=0.3, labelspacing=0.25, ncol=1)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"✓ wrote {OUT}  ({len(cats)} cell types, {adata.n_obs} cells)")


if __name__ == "__main__":
    main()
