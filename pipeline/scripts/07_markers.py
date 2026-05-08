#!/usr/bin/env python
"""Manual marker validation for cluster and annotation review."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


DEFAULT_MARKER_DB = "pipeline/config/cell_type_markers.tsv"
PREDICTION_COLUMNS = [
    "final_cell_type",
    "popv_prediction",
    "popv_majority_vote_prediction",
    "scanvi_label",
    "paul15_label",
]
NON_LINEAGE_MARKER_SETS = {"Low_Quality_Ambient", "Stress_Response", "Cycling"}


def normalize_label(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def read_marker_table(path: str) -> pd.DataFrame:
    markers = pd.read_csv(path, sep="\t")
    required = {"cell_type", "marker", "direction"}
    missing = required - set(markers.columns)
    if missing:
        raise ValueError(f"Marker table is missing columns: {sorted(missing)}")
    return markers


def marker_dict(adata: sc.AnnData, markers: pd.DataFrame, out: Path) -> dict[str, list[str]]:
    positive = markers.loc[markers["direction"].eq("positive")].copy()
    grouped = {}
    missing = []
    for cell_type, sub in positive.groupby("cell_type", sort=False):
        genes = [gene for gene in sub["marker"].astype(str) if gene in adata.var_names]
        if genes:
            grouped[cell_type] = genes
        missing.extend([gene for gene in sub["marker"].astype(str) if gene not in adata.var_names])
    if missing:
        pd.Series(sorted(set(missing)), name="missing_marker").to_csv(out / "missing_markers.csv", index=False)
        print(f"Missing markers: {sorted(set(missing))}")
    return grouped


def expression_frame(adata: sc.AnnData, genes: list[str], layer: str) -> pd.DataFrame:
    x = adata[:, genes].layers[layer] if layer in adata.layers else adata[:, genes].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return pd.DataFrame(x, index=adata.obs_names, columns=genes)


def score_marker_sets(adata: sc.AnnData, markers: dict[str, list[str]], groupby: str, layer: str) -> pd.DataFrame:
    rows = []
    groups = adata.obs[groupby].astype(str)
    for cell_type, genes in markers.items():
        expr = expression_frame(adata, genes, layer)
        per_cell_score = expr.mean(axis=1)
        scores = per_cell_score.groupby(groups, observed=True).mean()
        pct_expr = (expr.gt(0).mean(axis=1)).groupby(groups, observed=True).mean()
        for group, score in scores.items():
            rows.append(
                {
                    "groupby": groupby,
                    "group": str(group),
                    "marker_cell_type": cell_type,
                    "mean_marker_score": float(score),
                    "mean_fraction_markers_detected": float(pct_expr.loc[group]),
                    "n_markers_present": len(genes),
                    "markers_present": ",".join(genes),
                }
            )
    return pd.DataFrame(rows)


def plot_marker_heatmap(scores: pd.DataFrame, out: Path, filename: str) -> None:
    if scores.empty:
        return
    matrix = scores.pivot(index="group", columns="marker_cell_type", values="mean_marker_score").fillna(0)
    if matrix.shape[0] == 0 or matrix.shape[1] == 0:
        return
    height = max(4, min(18, 0.35 * matrix.shape[0] + 2))
    width = max(8, min(24, 0.35 * matrix.shape[1] + 4))
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(matrix, cmap="viridis", ax=ax)
    ax.set_xlabel("Marker set")
    ax.set_ylabel(scores["groupby"].iloc[0])
    plt.tight_layout()
    plt.savefig(out / filename, dpi=150)
    plt.close(fig)


def best_marker_matches(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame()
    rows = []
    for group, sub in scores.groupby("group", sort=False):
        top = sub.sort_values("mean_marker_score", ascending=False).head(5)
        for rank, (_, row) in enumerate(top.iterrows(), start=1):
            rows.append(
                {
                    "groupby": row["groupby"],
                    "group": group,
                    "rank": rank,
                    "marker_cell_type": row["marker_cell_type"],
                    "mean_marker_score": row["mean_marker_score"],
                    "mean_fraction_markers_detected": row["mean_fraction_markers_detected"],
                    "markers_present": row["markers_present"],
                }
            )
    return pd.DataFrame(rows)


def lineage_marker_scores(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return scores
    return scores.loc[~scores["marker_cell_type"].isin(NON_LINEAGE_MARKER_SETS)].copy()


def expected_marker_match(label: str, marker_cell_types: list[str]) -> str | None:
    label_norm = normalize_label(label)
    aliases = {
        "hsc": "HSC_MPP",
        "mpp": "HSC_MPP",
        "lmp": "MPP_LMPP",
        "lymph": "MPP_LMPP",
        "gmp": "CMP_GMP",
        "cmp": "CMP_GMP",
        "mep": "MEP_Erythroid",
        "ery": "MEP_Erythroid",
        "meg": "Megakaryocyte",
        "mk": "Megakaryocyte",
        "mega": "Megakaryocyte",
        "bas": "Basophil_Mast",
        "baso": "Basophil_Mast",
        "mast": "Basophil_Mast",
        "mono": "Monocyte",
        "dc": "DC",
        "neu": "Granulocyte_Neutrophil",
        "gran": "Granulocyte_Neutrophil",
        "bcell": "B_cell",
        "b": "B_cell",
        "tcell": "T_cell",
        "t": "T_cell",
        "nk": "NK_cell",
        "plasma": "Plasma_cell",
    }
    for key, marker_type in aliases.items():
        if key in label_norm and marker_type in marker_cell_types:
            return marker_type
    for marker_type in marker_cell_types:
        marker_norm = normalize_label(marker_type)
        if marker_norm in label_norm or label_norm in marker_norm:
            return marker_type
    return None


def write_prediction_marker_consistency(
    adata: sc.AnnData,
    out: Path,
    markers: dict[str, list[str]],
    layer: str,
) -> None:
    rows = []
    marker_types = list(markers.keys())
    for label_col in PREDICTION_COLUMNS:
        if label_col not in adata.obs:
            continue
        scores = score_marker_sets(adata, markers, label_col, layer)
        scores.to_csv(out / f"marker_scores_by_{label_col}.csv", index=False)
        plot_marker_heatmap(scores, out, f"heatmap_marker_scores_by_{label_col}.png")
        best = best_marker_matches(scores)
        best.to_csv(out / f"top_marker_matches_by_{label_col}.csv", index=False)
        best_lineage = best_marker_matches(lineage_marker_scores(scores))
        best_lineage.to_csv(out / f"top_lineage_marker_matches_by_{label_col}.csv", index=False)
        best_top = best_lineage.loc[best_lineage["rank"].eq(1)].set_index("group") if not best_lineage.empty else pd.DataFrame()
        for label in adata.obs[label_col].astype(str).unique():
            expected = expected_marker_match(label, marker_types)
            top_marker = None
            top_score = np.nan
            expected_score = np.nan
            if label in best_top.index:
                top_marker = best_top.loc[label, "marker_cell_type"]
                top_score = best_top.loc[label, "mean_marker_score"]
            if expected is not None:
                match = scores.loc[
                    scores["group"].eq(label) & scores["marker_cell_type"].eq(expected),
                    "mean_marker_score",
                ]
                if not match.empty:
                    expected_score = float(match.iloc[0])
            rows.append(
                {
                    "prediction_column": label_col,
                    "predicted_label": label,
                    "expected_marker_set": expected if expected is not None else "unmatched_label",
                    "top_marker_set": top_marker,
                    "top_marker_score": top_score,
                    "expected_marker_score": expected_score,
                    "expected_is_top_marker_set": expected is not None and expected == top_marker,
                    "n_cells": int(adata.obs[label_col].astype(str).eq(label).sum()),
                }
            )
    pd.DataFrame(rows).to_csv(out / "annotation_marker_consistency.csv", index=False)


def parse_watchlist(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def write_qc_summaries(adata: sc.AnnData, out: Path, cluster_key: str, watchlist: list[str]) -> None:
    qc_cols = [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
        "doubletfinder_score",
    ]
    qc_cols = [col for col in qc_cols if col in adata.obs]
    if not qc_cols or cluster_key not in adata.obs:
        return
    adata.obs.groupby(cluster_key, observed=True)[qc_cols].agg(["count", "median", "mean"]).to_csv(
        out / "cluster_qc_summary.csv"
    )
    clusters = adata.obs[cluster_key].astype(str)
    watch_mask = clusters.isin(watchlist)
    if watch_mask.any():
        adata.obs["watchlist_cluster"] = "other"
        adata.obs.loc[watch_mask, "watchlist_cluster"] = clusters.loc[watch_mask].map(lambda x: f"cluster_{x}")
        adata.obs.loc[watch_mask].groupby(cluster_key, observed=True)[qc_cols].agg(["count", "median", "mean"]).to_csv(
            out / "watchlist_qc_summary.csv"
        )
        sc.pl.umap(adata, color="watchlist_cluster", show=False)
        plt.savefig(out / "umap_watchlist_clusters.png", dpi=150, bbox_inches="tight")
        plt.close()


def write_cluster_composition(adata: sc.AnnData, out: Path, cluster_key: str) -> None:
    if cluster_key not in adata.obs:
        return
    clusters = adata.obs[cluster_key].astype(str)
    rows = []
    for variable in ["pool", "sample_id", "population", "donor", "treatment"]:
        if variable not in adata.obs:
            continue
        counts = pd.crosstab(clusters, adata.obs[variable].astype(str))
        proportions = counts.div(counts.sum(axis=1), axis=0).fillna(0)
        counts.to_csv(out / f"{cluster_key}_{variable}_counts.csv")
        proportions.to_csv(out / f"{cluster_key}_{variable}_proportions.csv")
        for cluster in counts.index:
            rows.append(
                {
                    "cluster": cluster,
                    "variable": variable,
                    "dominant_category": proportions.loc[cluster].idxmax(),
                    "dominant_fraction": proportions.loc[cluster].max(),
                }
            )
    pd.DataFrame(rows).to_csv(out / "cluster_composition_dominance.csv", index=False)


def main(in_path: str, out_dir: str, cluster_key: str, marker_db: str, layer: str, watchlist: list[str]) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)
    markers_table = read_marker_table(marker_db)
    markers = marker_dict(adata, markers_table, out)
    pd.DataFrame(
        [{"cell_type": cell_type, "markers_present": ",".join(genes), "n_markers_present": len(genes)}
         for cell_type, genes in markers.items()]
    ).to_csv(out / "marker_sets_present.csv", index=False)

    write_qc_summaries(adata, out, cluster_key, watchlist)
    write_cluster_composition(adata, out, cluster_key)

    if cluster_key in adata.obs:
        cluster_scores = score_marker_sets(adata, markers, cluster_key, layer)
        cluster_scores.to_csv(out / "marker_scores_by_cluster.csv", index=False)
        best_marker_matches(cluster_scores).to_csv(out / "top_marker_matches_by_cluster.csv", index=False)
        best_marker_matches(lineage_marker_scores(cluster_scores)).to_csv(
            out / "top_lineage_marker_matches_by_cluster.csv",
            index=False,
        )
        plot_marker_heatmap(cluster_scores, out, "heatmap_marker_scores_by_cluster.png")
        sc.pl.dotplot(
            adata,
            var_names=markers,
            groupby=cluster_key,
            use_raw=False,
            layer=layer if layer in adata.layers else None,
            standard_scale="var",
            show=False,
        )
        plt.savefig(out / "dotplot_clusters_markers.png", dpi=150, bbox_inches="tight")
        plt.close()

    # Per-prediction-column dotplots — intermediate ML outputs (popv, scanvi, etc.), not final annotations.
    # for label_col in PREDICTION_COLUMNS:
    #     if label_col in adata.obs:
    #         sc.pl.dotplot(
    #             adata,
    #             var_names=markers,
    #             groupby=label_col,
    #             use_raw=False,
    #             layer=layer if layer in adata.layers else None,
    #             standard_scale="var",
    #             show=False,
    #         )
    #         plt.savefig(out / f"dotplot_{label_col}_markers.png", dpi=150, bbox_inches="tight")
    #         plt.close()

    write_prediction_marker_consistency(adata, out, markers, layer)

    # Per-gene expression UMAPs — exploratory during annotation; generates 10 files per run.
    # for gene in ["Kit", "Ly6a", "Hlf", "Mecom", "Mki67", "Cd79a", "Cd3d", "Gata1", "Mpo", "Csf3r"]:
    #     if gene in adata.var_names:
    #         sc.pl.umap(adata, color=gene, layer=layer if layer in adata.layers else None, show=False, vmax="p99")
    #         plt.savefig(out / f"umap_expr_{gene}.png", dpi=150, bbox_inches="tight")
    #         plt.close()

    if cluster_key in adata.obs:
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon", layer=layer, use_raw=False)
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
        plt.savefig(out / "dotplot_top_de_genes.png", dpi=150, bbox_inches="tight")
        plt.close()
        de = sc.get.rank_genes_groups_df(adata, group=None)
        de.to_csv(out / "de_genes_per_cluster.csv", index=False)
        if watchlist:
            de.loc[de["group"].astype(str).isin(watchlist)].to_csv(out / "watchlist_de_genes.csv", index=False)
    print(f"Saved marker validation outputs to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cluster_key", default="leiden_r0.5")
    parser.add_argument("--marker-db", default=DEFAULT_MARKER_DB)
    parser.add_argument("--layer", default="lognorm")
    parser.add_argument("--watchlist-clusters", default="10,12")
    args = parser.parse_args()
    main(args.input, args.out, args.cluster_key, args.marker_db, args.layer, parse_watchlist(args.watchlist_clusters))
