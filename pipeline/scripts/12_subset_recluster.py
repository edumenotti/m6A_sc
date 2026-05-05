#!/usr/bin/env python
"""Subset reclustering for B cell and erythroid compartments.

Extracts leiden_r1.0 clusters 15 (B cell) and 13 (erythroid) from the
annotated object, re-embeds at higher resolution, scores lineage markers,
assigns level2 labels, and updates manual_level2 in the full object.
If all subclusters call the same label, keeps the broad label unchanged.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 1

B_CLUSTER = "15"
ERY_CLUSTER = "13"
CLUSTER_KEY = "leiden_r1.0"

PRO_PRE_B    = ["Vpreb1", "Vpreb3", "Ebf1", "Dntt", "Bach2"]
MATURE_B     = ["H2-Aa", "H2-Eb1", "H2-Ab1", "Bank1", "Ms4a1", "Cd74"]
MEP_ERY      = ["Gata1", "Klf1", "Itga2b", "Gata2", "Hba-a1"]
ERYTHROBLAST = ["Car2", "Car1", "Blvrb", "Aqp1", "Ermap", "Slc4a1"]


def present(adata: sc.AnnData, genes: list[str]) -> list[str]:
    return [g for g in genes if g in adata.var_names]


def mean_expr(adata_sub: sc.AnnData, genes: list[str]) -> float:
    ok = present(adata_sub, genes)
    if not ok:
        return 0.0
    arr = adata_sub[:, ok].X
    if hasattr(arr, "toarray"):
        arr = arr.toarray()
    return float(np.asarray(arr).mean())


def recluster(adata_sub: sc.AnnData, res: float = 0.5,
              n_hvg: int = 400, n_pcs: int = 15, n_neighbors: int = 10) -> sc.AnnData:
    """Re-embed: set X=lognorm → HVG → PCA → neighbors → leiden → UMAP."""
    if "lognorm" in adata_sub.layers:
        adata_sub.X = adata_sub.layers["lognorm"].copy()
    sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_hvg, subset=False)
    sc.pp.pca(adata_sub, n_comps=n_pcs, use_highly_variable=True)
    sc.pp.neighbors(adata_sub, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(adata_sub, resolution=res, key_added="leiden_subset")
    sc.tl.umap(adata_sub)
    return adata_sub


def make_calls(adata_sub: sc.AnnData, pos_genes: list[str], neg_genes: list[str],
               pos_label: str, neg_label: str) -> dict[str, str]:
    calls: dict[str, str] = {}
    for c in sorted(adata_sub.obs["leiden_subset"].unique(), key=int):
        sub = adata_sub[adata_sub.obs["leiden_subset"] == c]
        pos = mean_expr(sub, pos_genes)
        neg = mean_expr(sub, neg_genes)
        calls[c] = pos_label if pos >= neg else neg_label
    return calls


def plot_subset(adata_sub: sc.AnnData, cols: list[str], out: Path, prefix: str) -> None:
    for col in cols:
        in_obs = col in adata_sub.obs.columns
        in_var = col in adata_sub.var_names
        if not in_obs and not in_var:
            continue
        sc.pl.umap(adata_sub, color=col, show=False, frameon=False,
                   legend_loc="right margin", legend_fontsize=7)
        plt.savefig(out / f"{prefix}_umap_{col.replace('-', '_').replace('/', '_')}.png",
                    dpi=150, bbox_inches="tight")
        plt.close()


def process_compartment(adata: sc.AnnData, cluster_id: str,
                        pos_genes: list[str], neg_genes: list[str],
                        pos_label: str, neg_label: str, broad_label: str,
                        plot_genes: list[str], prefix: str, out: Path) -> pd.Series:
    mask = adata.obs[CLUSTER_KEY].astype(str) == cluster_id
    adata_sub = adata[mask].copy()
    print(f"  {prefix} subset: {adata_sub.n_obs} cells")

    adata_sub = recluster(adata_sub)
    calls = make_calls(adata_sub, pos_genes, neg_genes, pos_label, neg_label)

    rows = []
    for c, label in calls.items():
        sub = adata_sub[adata_sub.obs["leiden_subset"] == c]
        rows.append({
            "compartment": prefix,
            "subcluster": c,
            "n_cells": int(sub.n_obs),
            "pos_markers_mean": round(mean_expr(sub, pos_genes), 4),
            "neg_markers_mean": round(mean_expr(sub, neg_genes), 4),
            "level2_call": label,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out / f"{prefix}_subcluster_summary.csv", index=False)
    print(summary[["subcluster", "n_cells", "pos_markers_mean", "neg_markers_mean", "level2_call"]].to_string(index=False))

    if summary["level2_call"].nunique() == 1:
        print(f"  [{prefix}] No split supported — keeping {broad_label}")
        result = pd.Series(broad_label, index=adata_sub.obs_names, name="level2_new")
    else:
        result = adata_sub.obs["leiden_subset"].map(calls).fillna(broad_label)
        result.name = "level2_new"
        counts = summary.groupby("level2_call")["n_cells"].sum().to_dict()
        print(f"  [{prefix}] Split: {counts}")

    plot_subset(adata_sub, ["leiden_subset", "level2_new"] + plot_genes, out, prefix)
    return result


def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells")

    b_updates = process_compartment(
        adata, B_CLUSTER,
        pos_genes=PRO_PRE_B, neg_genes=MATURE_B,
        pos_label="Pro_Pre_B", neg_label="Mature_B", broad_label="B_cell_broad",
        plot_genes=["Vpreb3", "H2-Aa"],
        prefix="b_cell", out=out,
    )

    ery_updates = process_compartment(
        adata, ERY_CLUSTER,
        pos_genes=MEP_ERY, neg_genes=ERYTHROBLAST,
        pos_label="MEP_Erythroid", neg_label="Erythroblast", broad_label="Erythroid_broad",
        plot_genes=["Gata1", "Car2"],
        prefix="erythroid", out=out,
    )

    all_updates = pd.concat([b_updates, ery_updates])
    adata.obs["manual_level2"] = adata.obs["manual_level2"].astype(str)
    adata.obs.loc[all_updates.index, "manual_level2"] = all_updates.values
    adata.obs["manual_level2"] = adata.obs["manual_level2"].astype("category")

    counts = (
        adata.obs.groupby(["manual_level1", "manual_level2"], observed=True)
        .size().reset_index(name="n_cells")
    )
    counts["pct"] = (counts["n_cells"] / adata.n_obs * 100).round(2)
    counts.to_csv(out / "final_level1_level2_counts.csv", index=False)

    adata.write_h5ad(out / "adata_annotated_final.h5ad")
    print(f"Saved {out}/adata_annotated_final.h5ad")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/11_annotation_r1.0/adata_manual_level1.h5ad")
    parser.add_argument("--out",   default="results/12_subset_recluster")
    main(parser.parse_args())
