#!/usr/bin/env python
"""Score manual marker sets and generate review plots before final annotation."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


DEFAULT_MARKER_DB = "pipeline/config/manual_annotation_markers_skull_immune.tsv"
DEFAULT_INPUT = "results/08_reconcile_annotations/adata_reconciled.h5ad"
QC_LEVEL1 = {"QC"}
LOW_SPECIFICITY_LEVEL2 = {"Pre_neutrophil"}


def read_markers(path: str) -> pd.DataFrame:
    markers = pd.read_csv(path, sep="\t")
    required = {"level1", "level2", "marker", "role", "source", "priority", "note"}
    missing = required - set(markers.columns)
    if missing:
        raise ValueError(f"Marker table is missing columns: {sorted(missing)}")
    markers = markers.loc[markers["role"].eq("positive")].copy()
    for col in ["level1", "level2", "marker", "priority"]:
        markers[col] = markers[col].astype(str)
    return markers


def present_marker_table(adata: sc.AnnData, markers: pd.DataFrame, out: Path) -> pd.DataFrame:
    present = markers.loc[markers["marker"].isin(adata.var_names)].copy()
    missing = markers.loc[~markers["marker"].isin(adata.var_names)].copy()
    present.to_csv(out / "manual_markers_present.tsv", sep="\t", index=False)
    missing.to_csv(out / "manual_markers_missing.tsv", sep="\t", index=False)
    return present


def marker_sets(markers: pd.DataFrame, group_col: str, max_markers: int | None = None) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = {}
    priority_order = {"high": 0, "medium": 1, "low": 2}
    ordered = markers.assign(_priority=markers["priority"].map(priority_order).fillna(9))
    ordered = ordered.sort_values([group_col, "_priority", "marker"])
    for group, sub in ordered.groupby(group_col, sort=False):
        genes = list(dict.fromkeys(sub["marker"].astype(str)))
        if max_markers is not None:
            genes = genes[:max_markers]
        if genes:
            grouped[str(group)] = genes
    return grouped


def expression_frame(adata: sc.AnnData, genes: list[str], layer: str) -> pd.DataFrame:
    x = adata[:, genes].layers[layer] if layer in adata.layers else adata[:, genes].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return pd.DataFrame(x, index=adata.obs_names, columns=genes)


def score_sets(adata: sc.AnnData, sets: dict[str, list[str]], groupby: str, layer: str, set_kind: str) -> pd.DataFrame:
    groups = adata.obs[groupby].astype(str)
    rows = []
    for marker_set, genes in sets.items():
        expr = expression_frame(adata, genes, layer)
        per_cell_score = expr.mean(axis=1)
        per_cell_fraction = expr.gt(0).mean(axis=1)
        mean_scores = per_cell_score.groupby(groups, observed=True).mean()
        mean_fraction = per_cell_fraction.groupby(groups, observed=True).mean()
        for group, score in mean_scores.items():
            rows.append(
                {
                    "groupby": groupby,
                    "group": str(group),
                    "set_kind": set_kind,
                    "marker_set": marker_set,
                    "mean_marker_score": float(score),
                    "mean_fraction_markers_detected": float(mean_fraction.loc[group]),
                    "n_markers_present": len(genes),
                    "markers_present": ",".join(genes),
                }
            )
    return pd.DataFrame(rows)


def best_matches(scores: pd.DataFrame, exclude_sets: set[str] | None = None, top_n: int = 5) -> pd.DataFrame:
    if scores.empty:
        return scores
    work = scores.copy()
    if exclude_sets:
        work = work.loc[~work["marker_set"].isin(exclude_sets)].copy()
    rows = []
    for group, sub in work.groupby("group", sort=False):
        top = sub.sort_values("mean_marker_score", ascending=False).head(top_n)
        for rank, (_, row) in enumerate(top.iterrows(), start=1):
            rows.append(
                {
                    "groupby": row["groupby"],
                    "group": group,
                    "set_kind": row["set_kind"],
                    "rank": rank,
                    "marker_set": row["marker_set"],
                    "mean_marker_score": row["mean_marker_score"],
                    "mean_fraction_markers_detected": row["mean_fraction_markers_detected"],
                    "n_markers_present": row["n_markers_present"],
                    "markers_present": row["markers_present"],
                }
            )
    return pd.DataFrame(rows)


def best_level2_within_level1(
    level2_scores: pd.DataFrame,
    level1_best: pd.DataFrame,
    level2_to_level1: dict[str, str],
) -> pd.DataFrame:
    if level2_scores.empty or level1_best.empty:
        return pd.DataFrame()
    top_level1 = level1_best.loc[level1_best["rank"].eq(1)].set_index("group")["marker_set"].to_dict()
    rows = []
    for group, sub in level2_scores.groupby("group", sort=False):
        dominant_level1 = top_level1.get(group)
        if dominant_level1 is None:
            continue
        candidates = sub.loc[sub["marker_set"].map(level2_to_level1).eq(dominant_level1)].copy()
        if candidates.empty:
            continue
        strict = candidates.loc[~candidates["marker_set"].isin(LOW_SPECIFICITY_LEVEL2)].copy()
        ranking_pool = strict if not strict.empty else candidates
        top = ranking_pool.sort_values("mean_marker_score", ascending=False).head(5)
        for rank, (_, row) in enumerate(top.iterrows(), start=1):
            rows.append(
                {
                    "groupby": row["groupby"],
                    "group": group,
                    "set_kind": row["set_kind"],
                    "rank": rank,
                    "dominant_level1": dominant_level1,
                    "marker_set": row["marker_set"],
                    "mean_marker_score": row["mean_marker_score"],
                    "mean_fraction_markers_detected": row["mean_fraction_markers_detected"],
                    "n_markers_present": row["n_markers_present"],
                    "markers_present": row["markers_present"],
                    "selection_note": "low_specificity_level2_excluded" if not strict.empty else "only_low_specificity_candidates",
                }
            )
    return pd.DataFrame(rows)


def plot_heatmap(scores: pd.DataFrame, out_file: Path, title: str) -> None:
    if scores.empty:
        return
    matrix = scores.pivot(index="group", columns="marker_set", values="mean_marker_score").fillna(0)
    try:
        matrix = matrix.loc[sorted(matrix.index, key=lambda x: int(x) if str(x).isdigit() else str(x))]
    except TypeError:
        pass
    height = max(6, min(18, 0.35 * matrix.shape[0] + 2))
    width = max(10, min(36, 0.32 * matrix.shape[1] + 4))
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(matrix, cmap="viridis", ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Marker set")
    ax.set_ylabel(scores["groupby"].iloc[0])
    plt.tight_layout()
    fig.savefig(out_file, dpi=150)
    plt.close(fig)


def plot_dotplot(
    adata: sc.AnnData,
    var_names: dict[str, list[str]],
    groupby: str,
    out_file: Path,
    layer: str,
    title: str,
) -> None:
    if not var_names:
        return
    width = max(8, min(44, sum(len(v) for v in var_names.values()) * 0.32 + 4))
    height = max(5, min(20, adata.obs[groupby].astype(str).nunique() * 0.35 + 3))
    sc.pl.dotplot(
        adata,
        var_names=var_names,
        groupby=groupby,
        use_raw=False,
        layer=layer if layer in adata.layers else None,
        standard_scale="var",
        show=False,
    )
    fig = plt.gcf()
    fig.set_size_inches(width, height)
    fig.suptitle(title, y=1.02)
    fig.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close(fig)


def top_categories(series: pd.Series, n: int = 5) -> str:
    counts = series.astype(str).value_counts().head(n)
    return ";".join(f"{idx}:{val}" for idx, val in counts.items())


def de_top_genes(adata: sc.AnnData, cluster_key: str, layer: str, out: Path) -> dict[str, str]:
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon", layer=layer, use_raw=False)
    de = sc.get.rank_genes_groups_df(adata, group=None)
    de.to_csv(out / "de_genes_per_cluster.csv", index=False)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
    plt.savefig(out / "dotplot_top_de_genes_by_cluster.png", dpi=150, bbox_inches="tight")
    plt.close()
    return de.groupby("group")["names"].apply(lambda s: ",".join(s.head(8))).to_dict()


def write_cluster_review(
    adata: sc.AnnData,
    cluster_key: str,
    level1_best: pd.DataFrame,
    level2_best: pd.DataFrame,
    level2_within_level1_best: pd.DataFrame,
    top_de: dict[str, str],
    out: Path,
) -> pd.DataFrame:
    level1_top = level1_best.loc[level1_best["rank"].eq(1)].set_index("group")
    level2_top = level2_best.loc[level2_best["rank"].eq(1)].set_index("group")
    level2_within_top = level2_within_level1_best.loc[level2_within_level1_best["rank"].eq(1)].set_index("group")
    qc_cols = [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
        "doubletfinder_score",
    ]
    qc_cols = [col for col in qc_cols if col in adata.obs]
    rows = []
    obs = adata.obs.copy()
    obs[cluster_key] = obs[cluster_key].astype(str)
    for cluster, sub in obs.groupby(cluster_key, observed=True):
        row = {
            "cluster_key": cluster_key,
            "cluster": str(cluster),
            "n_cells": int(sub.shape[0]),
            "top_manual_level1": level1_top.loc[str(cluster), "marker_set"] if str(cluster) in level1_top.index else "",
            "top_manual_level1_score": level1_top.loc[str(cluster), "mean_marker_score"] if str(cluster) in level1_top.index else np.nan,
            "top_manual_level1_fraction_detected": level1_top.loc[str(cluster), "mean_fraction_markers_detected"] if str(cluster) in level1_top.index else np.nan,
            "top_manual_level2": level2_top.loc[str(cluster), "marker_set"] if str(cluster) in level2_top.index else "",
            "top_manual_level2_score": level2_top.loc[str(cluster), "mean_marker_score"] if str(cluster) in level2_top.index else np.nan,
            "top_manual_level2_fraction_detected": level2_top.loc[str(cluster), "mean_fraction_markers_detected"] if str(cluster) in level2_top.index else np.nan,
            "top_manual_level2_markers": level2_top.loc[str(cluster), "markers_present"] if str(cluster) in level2_top.index else "",
            "top_level2_within_level1": level2_within_top.loc[str(cluster), "marker_set"] if str(cluster) in level2_within_top.index else "",
            "top_level2_within_level1_score": level2_within_top.loc[str(cluster), "mean_marker_score"] if str(cluster) in level2_within_top.index else np.nan,
            "top_level2_within_level1_fraction_detected": level2_within_top.loc[str(cluster), "mean_fraction_markers_detected"] if str(cluster) in level2_within_top.index else np.nan,
            "top_level2_within_level1_markers": level2_within_top.loc[str(cluster), "markers_present"] if str(cluster) in level2_within_top.index else "",
            "level2_selection_note": level2_within_top.loc[str(cluster), "selection_note"] if str(cluster) in level2_within_top.index else "",
            "top_de_genes": top_de.get(str(cluster), ""),
        }
        for col in ["leiden_r0.5", "final_cell_type", "popv_prediction", "scanvi_label", "sample_id", "population"]:
            if col in sub:
                row[f"top_{col}_counts"] = top_categories(sub[col])
        for col in qc_cols:
            row[f"median_{col}"] = float(pd.to_numeric(sub[col], errors="coerce").median())
        rows.append(row)
    review = pd.DataFrame(rows)
    review["_sort"] = review["cluster"].map(lambda x: int(x) if str(x).isdigit() else 10**9)
    review = review.sort_values(["_sort", "cluster"]).drop(columns="_sort")
    review.to_csv(out / f"manual_cluster_review_{cluster_key}.csv", index=False)
    return review


def write_checkpoint(review: pd.DataFrame, out: Path, cluster_key: str) -> None:
    lines = [
        "# Manual Annotation Checkpoint",
        "",
        f"Cluster key: `{cluster_key}`",
        "",
        "This checkpoint is for human review before applying manual labels.",
        "",
        "## Top Marker Calls",
        "",
        "| cluster | n_cells | top level1 | top level2 | top DE genes | notes |",
        "|---|---:|---|---|---|---|",
    ]
    for _, row in review.iterrows():
        de = str(row["top_de_genes"]).replace("|", "/")
        lines.append(
            f"| {row['cluster']} | {row['n_cells']} | {row['top_manual_level1']} "
            f"({row['top_manual_level1_score']:.2f}) | {row['top_level2_within_level1']} "
            f"({row['top_level2_within_level1_score']:.2f}) | {de} | {row['level2_selection_note']} |"
        )
    lines.extend(
        [
            "",
            "## Files To Inspect",
            "",
            "- `dotplot_level1_core_markers_by_cluster.png`",
            "- `dotplot_level2_core_markers_by_cluster.png`",
            "- `dotplot_<level1>_level2_markers_by_cluster.png` for compartment-specific checks",
            f"- `manual_cluster_review_{cluster_key}.csv`",
            "- `top_level1_marker_matches_by_cluster.csv`",
            "- `top_level2_marker_matches_by_cluster.csv`",
            "- `top_level2_within_level1_marker_matches_by_cluster.csv`",
            "",
            "Use this file to record corrections before creating the manual mapping TSV.",
        ]
    )
    (out / "CHECKPOINT_manual_annotation_review.md").write_text("\n".join(lines) + "\n")


def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.input)
    markers = read_markers(args.marker_db)
    present = present_marker_table(adata, markers, out)

    level1_sets = marker_sets(present.loc[~present["level1"].isin(QC_LEVEL1)], "level1")
    level2_sets = marker_sets(present, "level2")
    level2_to_level1 = present.drop_duplicates("level2").set_index("level2")["level1"].to_dict()
    level1_core = marker_sets(
        present.loc[~present["level1"].isin(QC_LEVEL1) & present["priority"].isin(["high", "medium"])],
        "level1",
        max_markers=args.max_markers_per_level1,
    )
    level2_core = marker_sets(
        present.loc[present["priority"].isin(["high", "medium"])],
        "level2",
        max_markers=args.max_markers_per_level2,
    )

    pd.DataFrame(
        [{"set_kind": "level1", "marker_set": k, "markers_present": ",".join(v), "n_markers_present": len(v)} for k, v in level1_sets.items()]
        + [{"set_kind": "level2", "marker_set": k, "markers_present": ",".join(v), "n_markers_present": len(v)} for k, v in level2_sets.items()]
    ).to_csv(out / "manual_marker_sets_present.csv", index=False)

    level1_scores = score_sets(adata, level1_sets, args.cluster_key, args.layer, "level1")
    level2_scores = score_sets(adata, level2_sets, args.cluster_key, args.layer, "level2")
    level1_scores.to_csv(out / "manual_level1_marker_scores_by_cluster.csv", index=False)
    level2_scores.to_csv(out / "manual_level2_marker_scores_by_cluster.csv", index=False)
    level1_best = best_matches(level1_scores)
    level2_best = best_matches(level2_scores, exclude_sets={"Cycling", "Stress_Response", "Low_Quality_Ambient"})
    level2_within_level1_best = best_level2_within_level1(level2_scores, level1_best, level2_to_level1)
    level1_best.to_csv(out / "top_level1_marker_matches_by_cluster.csv", index=False)
    level2_best.to_csv(out / "top_level2_marker_matches_by_cluster.csv", index=False)
    level2_within_level1_best.to_csv(out / "top_level2_within_level1_marker_matches_by_cluster.csv", index=False)

    plot_heatmap(level1_scores, out / "heatmap_level1_marker_scores_by_cluster.png", "Manual level1 marker scores")
    plot_heatmap(level2_scores, out / "heatmap_level2_marker_scores_by_cluster.png", "Manual level2 marker scores")
    plot_dotplot(
        adata,
        level1_core,
        args.cluster_key,
        out / "dotplot_level1_core_markers_by_cluster.png",
        args.layer,
        "Manual level1 core markers",
    )
    plot_dotplot(
        adata,
        level2_core,
        args.cluster_key,
        out / "dotplot_level2_core_markers_by_cluster.png",
        args.layer,
        "Manual level2 core markers",
    )

    for level1, sub in present.groupby("level1", sort=False):
        sets = marker_sets(sub.loc[sub["priority"].isin(["high", "medium"])], "level2", max_markers=args.max_markers_per_level2)
        if sets:
            safe = str(level1).replace("/", "_").replace(" ", "_")
            plot_dotplot(
                adata,
                sets,
                args.cluster_key,
                out / f"dotplot_{safe}_level2_markers_by_cluster.png",
                args.layer,
                f"{level1} level2 markers",
            )

    # Per-gene expression UMAPs — exploratory during manual review; generates ~23 files per run.
    # for gene in args.umap_genes.split(","):
    #     gene = gene.strip()
    #     if gene and gene in adata.var_names:
    #         sc.pl.umap(adata, color=gene, layer=args.layer if args.layer in adata.layers else None, show=False, vmax="p99")
    #         plt.savefig(out / f"umap_expr_{gene}.png", dpi=150, bbox_inches="tight")
    #         plt.close()

    top_de = de_top_genes(adata, args.cluster_key, args.layer, out)
    review = write_cluster_review(adata, args.cluster_key, level1_best, level2_best, level2_within_level1_best, top_de, out)
    write_checkpoint(review, out, args.cluster_key)
    print(f"Saved manual marker review outputs to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--out", required=True)
    parser.add_argument("--marker-db", default=DEFAULT_MARKER_DB)
    parser.add_argument("--cluster-key", default="leiden_r1.0")
    parser.add_argument("--layer", default="lognorm")
    parser.add_argument("--max-markers-per-level1", type=int, default=10)
    parser.add_argument("--max-markers-per-level2", type=int, default=4)
    parser.add_argument(
        "--umap-genes",
        default="Hlf,Kit,Ly6a,Cd34,Mpo,Csf3r,Ltf,Ngp,Retnlg,Cxcr2,Ccr2,Ly6c2,Vcan,Cd79a,Ighm,Cd3d,Jchain,Xbp1,Gata1,Pf4,Flt3,Irf8,Siglech",
    )
    main(parser.parse_args())
