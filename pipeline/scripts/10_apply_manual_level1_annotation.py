#!/usr/bin/env python
"""Apply an explicit manual level1 annotation map to an AnnData object."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


DEFAULT_INPUT = "results/08_reconcile_annotations/adata_reconciled.h5ad"
DEFAULT_MAP = "pipeline/config/manual_annotation_level1_map_leiden_r2.0.tsv"


def read_annotation_map(path: str, cluster_key: str) -> pd.DataFrame:
    mapping = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {
        "cluster_key",
        "cluster",
        "manual_level1",
        "manual_level2",
        "manual_annotation_confidence",
        "manual_review_flag",
        "manual_annotation_basis",
    }
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(f"Manual map is missing columns: {sorted(missing)}")

    mapping = mapping.loc[mapping["cluster_key"].eq(cluster_key)].copy()
    if mapping.empty:
        raise ValueError(f"No rows for cluster_key={cluster_key!r} in {path}")
    duplicated = mapping["cluster"].duplicated(keep=False)
    if duplicated.any():
        dupes = sorted(mapping.loc[duplicated, "cluster"].unique())
        raise ValueError(f"Duplicate manual map rows for clusters: {dupes}")
    return mapping


def top_categories(series: pd.Series, n: int = 5) -> str:
    counts = series.astype(str).value_counts().head(n)
    return ";".join(f"{idx}:{val}" for idx, val in counts.items())


def validate_clusters(adata: sc.AnnData, mapping: pd.DataFrame, cluster_key: str) -> None:
    if cluster_key not in adata.obs:
        raise ValueError(f"{cluster_key!r} is not present in adata.obs")
    observed = set(adata.obs[cluster_key].astype(str).unique())
    mapped = set(mapping["cluster"].astype(str))
    missing = sorted(observed - mapped, key=lambda x: int(x) if x.isdigit() else x)
    extra = sorted(mapped - observed, key=lambda x: int(x) if x.isdigit() else x)
    if missing:
        raise ValueError(f"Clusters missing from manual map: {missing}")
    if extra:
        raise ValueError(f"Manual map contains clusters not present in object: {extra}")


def apply_map(
    adata: sc.AnnData,
    mapping: pd.DataFrame,
    cluster_key: str,
    annotation_version: str,
) -> sc.AnnData:
    validate_clusters(adata, mapping, cluster_key)
    keyed = mapping.set_index("cluster")
    clusters = adata.obs[cluster_key].astype(str)

    adata.obs["manual_annotation_cluster_key"] = cluster_key
    adata.obs["manual_annotation_cluster"] = clusters.to_numpy()
    adata.obs["manual_level1"] = clusters.map(keyed["manual_level1"]).astype("category").to_numpy()
    adata.obs["manual_level2"] = clusters.map(keyed["manual_level2"]).astype("category").to_numpy()
    adata.obs["manual_annotation_confidence"] = clusters.map(keyed["manual_annotation_confidence"]).astype("category").to_numpy()
    adata.obs["manual_review_flag"] = clusters.map(keyed["manual_review_flag"]).astype("category").to_numpy()
    adata.obs["manual_annotation_basis"] = clusters.map(keyed["manual_annotation_basis"]).to_numpy()
    adata.obs["manual_annotation_source"] = "manual_marker_review"
    adata.obs["manual_annotation_version"] = annotation_version
    return adata


def write_cluster_summary(adata: sc.AnnData, mapping: pd.DataFrame, out: Path, cluster_key: str) -> None:
    obs = adata.obs.copy()
    obs[cluster_key] = obs[cluster_key].astype(str)
    rows = []
    for cluster, sub in obs.groupby(cluster_key, observed=True):
        row = {
            "cluster_key": cluster_key,
            "cluster": str(cluster),
            "n_cells": int(sub.shape[0]),
            "manual_level1": sub["manual_level1"].astype(str).iloc[0],
            "manual_level2": sub["manual_level2"].astype(str).iloc[0],
            "manual_annotation_confidence": sub["manual_annotation_confidence"].astype(str).iloc[0],
            "manual_review_flag": sub["manual_review_flag"].astype(str).iloc[0],
            "manual_annotation_basis": sub["manual_annotation_basis"].astype(str).iloc[0],
        }
        for col in ["final_cell_type", "popv_prediction", "scanvi_label", "population", "sample_id"]:
            if col in sub:
                row[f"top_{col}_counts"] = top_categories(sub[col])
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary["_sort"] = summary["cluster"].map(lambda x: int(x) if str(x).isdigit() else 10**9)
    summary = summary.sort_values(["_sort", "cluster"]).drop(columns="_sort")
    summary.to_csv(out / f"manual_level1_cluster_summary_{cluster_key}.csv", index=False)

    level1 = (
        obs["manual_level1"]
        .astype(str)
        .value_counts()
        .rename_axis("manual_level1")
        .reset_index(name="n_cells")
    )
    level1["fraction_cells"] = level1["n_cells"] / int(adata.n_obs)
    level1.to_csv(out / "manual_level1_counts.csv", index=False)

    if "final_cell_type" in obs:
        pd.crosstab(obs["manual_level1"].astype(str), obs["final_cell_type"].astype(str)).to_csv(
            out / "manual_level1_vs_reconciled_final_cell_type.csv"
        )
    if "population" in obs:
        pd.crosstab(obs["manual_level1"].astype(str), obs["population"].astype(str)).to_csv(
            out / "manual_level1_by_population.csv"
        )
    if "sample_id" in obs:
        pd.crosstab(obs["manual_level1"].astype(str), obs["sample_id"].astype(str)).to_csv(
            out / "manual_level1_by_sample_id.csv"
        )


def plot_umaps(adata: sc.AnnData, out: Path) -> None:
    if "X_umap" not in adata.obsm:
        return
    for col in ["manual_level1", "manual_level2", "manual_annotation_confidence", "manual_review_flag"]:
        sc.pl.umap(
            adata,
            color=col,
            show=False,
            legend_loc="right margin",
            frameon=False,
        )
        plt.savefig(out / f"umap_{col}.png", dpi=150, bbox_inches="tight")
        plt.close()


def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    mapping = read_annotation_map(args.map, args.cluster_key)
    adata = apply_map(adata, mapping, args.cluster_key, args.annotation_version)

    mapping.to_csv(out / Path(args.map).name, sep="\t", index=False)
    adata.obs[
        [
            args.cluster_key,
            "manual_level1",
            "manual_level2",
            "manual_annotation_confidence",
            "manual_review_flag",
            "manual_annotation_basis",
            "manual_annotation_source",
            "manual_annotation_version",
        ]
    ].to_csv(out / "manual_level1_per_cell.csv")
    write_cluster_summary(adata, mapping, out, args.cluster_key)
    plot_umaps(adata, out)
    adata.write_h5ad(out / "adata_manual_level1.h5ad")
    print(f"Saved manual level1 annotation outputs to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--map", default=DEFAULT_MAP)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cluster-key", default="leiden_r2.0")
    parser.add_argument("--annotation-version", default="manual_level1_r2.0_2026-05-05")
    main(parser.parse_args())
