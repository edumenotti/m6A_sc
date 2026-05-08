#!/usr/bin/env python
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


DEFAULT_RESOLUTIONS = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0]
COMPOSITION_KEYS = ["pool", "sample_id", "population", "donor", "treatment"]


def format_resolution(resolution: float) -> str:
    return str(float(resolution))


def parse_resolutions(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def write_cluster_composition(
    adata: sc.AnnData,
    out: Path,
    cluster_keys: list[str],
    composition_keys: list[str],
) -> None:
    summary_rows = []
    long_rows = []

    available_keys = [key for key in composition_keys if key in adata.obs]
    missing_keys = sorted(set(composition_keys) - set(available_keys))
    if missing_keys:
        print(f"Skipping missing composition keys: {', '.join(missing_keys)}")

    for cluster_key in cluster_keys:
        clusters = adata.obs[cluster_key].astype(str)
        cluster_sizes = clusters.value_counts().sort_index()

        summary_by_cluster = {
            cluster: {
                "cluster_key": cluster_key,
                "cluster": cluster,
                "n_cells": int(n_cells),
                "pct_cells": float(n_cells / adata.n_obs),
            }
            for cluster, n_cells in cluster_sizes.items()
        }

        for variable in available_keys:
            values = adata.obs[variable].astype(str)
            counts = pd.crosstab(clusters, values)
            counts.index.name = "cluster"
            counts.columns.name = variable
            counts.to_csv(out / f"{cluster_key}_{variable}_counts.csv")

            proportions = counts.div(counts.sum(axis=1), axis=0).fillna(0.0)
            proportions.to_csv(out / f"{cluster_key}_{variable}_proportions.csv")

            for cluster in counts.index:
                dominant_category = str(proportions.loc[cluster].idxmax())
                dominant_fraction = float(proportions.loc[cluster].max())
                summary_by_cluster[str(cluster)][f"dominant_{variable}"] = dominant_category
                summary_by_cluster[str(cluster)][f"dominant_{variable}_fraction"] = dominant_fraction

                for category, n_cells in counts.loc[cluster].items():
                    long_rows.append(
                        {
                            "cluster_key": cluster_key,
                            "cluster": str(cluster),
                            "variable": variable,
                            "category": str(category),
                            "n_cells": int(n_cells),
                            "prop_in_cluster": float(proportions.loc[cluster, category]),
                        }
                    )

        summary_rows.extend(summary_by_cluster.values())

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out / "cluster_dominance_summary.csv", index=False)
    pd.DataFrame(long_rows).to_csv(out / "cluster_composition_long.csv", index=False)

    resolution_rows = []
    for cluster_key, group in summary.groupby("cluster_key", sort=False):
        row = {
            "cluster_key": cluster_key,
            "n_clusters": int(group.shape[0]),
            "min_cluster_size": int(group["n_cells"].min()),
            "median_cluster_size": float(group["n_cells"].median()),
            "max_cluster_size": int(group["n_cells"].max()),
            "clusters_lt_100_cells": int((group["n_cells"] < 100).sum()),
        }
        for variable in available_keys:
            fraction_col = f"dominant_{variable}_fraction"
            row[f"median_dominant_{variable}_fraction"] = float(group[fraction_col].median())
            row[f"clusters_{variable}_gt80pct"] = int((group[fraction_col] >= 0.8).sum())
            row[f"clusters_{variable}_gt90pct"] = int((group[fraction_col] >= 0.9).sum())
        resolution_rows.append(row)

    pd.DataFrame(resolution_rows).to_csv(out / "cluster_resolution_summary.csv", index=False)


def main(in_path: str, out_dir: str, resolution: float, resolutions: list[float]) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)

    if "neighbors" not in adata.uns:
        if "X_scVI" not in adata.obsm:
            raise ValueError("Input AnnData has no neighbors graph and no X_scVI representation.")
        sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30, metric="euclidean")

    cluster_keys = []
    for res in resolutions:
        res_label = format_resolution(res)
        key = f"leiden_r{res_label}"
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=777)
        cluster_keys.append(key)
        sc.pl.umap(adata, color=key, legend_loc="on data", show=False)
        plt.savefig(out / f"umap_{key}.png", dpi=150, bbox_inches="tight")
        plt.close()

    # Batch/metadata UMAP panel — redundant with per-metadata UMAPs saved in 04_integrate.
    # metadata_colors = [key for key in ["population", "treatment", "donor", "pool"] if key in adata.obs]
    # sc.pl.umap(adata, color=metadata_colors, show=False)
    # plt.savefig(out / "umap_metadata.png", dpi=150, bbox_inches="tight")
    # plt.close()

    write_cluster_composition(adata, out, cluster_keys, COMPOSITION_KEYS)

    adata.write_h5ad(out / "adata_clustered.h5ad")
    key = f"leiden_r{format_resolution(resolution)}"
    print(f"Leiden r={resolution}: {adata.obs[key].nunique()} clusters")
    print(f"Saved: {out / 'adata_clustered.h5ad'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument(
        "--resolutions",
        default=",".join(str(res) for res in DEFAULT_RESOLUTIONS),
        help="Comma-separated Leiden resolutions to run.",
    )
    args = parser.parse_args()
    main(args.input, args.out, args.resolution, parse_resolutions(args.resolutions))
