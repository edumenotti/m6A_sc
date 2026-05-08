#!/usr/bin/env python
import argparse
from itertools import product
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


SAMPLE_MAP = {
    "1": ("D1", "DMSO", "LSK", "2"), "2": ("D1", "DMSO", "LSK", "1"),
    "3": ("D1", "DMSO", "LK", "2"), "4": ("D1", "DMSO", "LK", "1"),
    "5": ("D1", "STM", "LSK", "2"), "6": ("D1", "STM", "LSK", "1"),
    "7": ("D1", "STM", "LK", "2"), "8": ("D1", "STM", "LK", "1"),
    "9": ("D1", "STM", "I", "2"), "10": ("D1", "STM", "I", "1"),
    "11": ("D1", "DMSO", "I", "2"), "12": ("D1", "DMSO", "I", "1"),
    "13": ("D2", "STM", "LSK", "2"), "14": ("D2", "STM", "LSK", "1"),
    "15": ("D2", "DMSO", "LSK", "2"), "16": ("D2", "DMSO", "LSK", "1"),
    "17": ("D2", "STM", "LK", "2"), "18": ("D2", "STM", "LK", "1"),
    "19": ("D2", "DMSO", "LK", "2"), "20": ("D2", "DMSO", "LK", "1"),
    "21": ("D2", "STM", "I", "2"), "22": ("D2", "STM", "I", "1"),
    "23": ("D2", "DMSO", "I", "2"), "24": ("D2", "DMSO", "I", "1"),
}


def mad_bounds(values: np.ndarray, nmads: float, direction: str) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - med))
    if mad == 0 or np.isnan(mad):
        return -np.inf, np.inf
    low = med - nmads * mad if direction in ("lower", "both") else -np.inf
    high = med + nmads * mad if direction in ("upper", "both") else np.inf
    return low, high


def annotate_samples(adata: ad.AnnData) -> None:
    suffixes = [barcode.rsplit("-", 1)[-1] for barcode in adata.obs_names]
    meta = pd.DataFrame(
        [SAMPLE_MAP.get(suffix, ("unknown",) * 4) for suffix in suffixes],
        columns=["donor", "treatment", "population", "replicate"],
        index=adata.obs_names,
    )
    meta["sample_id"] = [
        f"{row.donor}_{row.treatment}_{row.population}_45_{row.replicate}"
        for row in meta.itertuples()
    ]
    adata.obs = meta


def load_qc(h5_path: str) -> ad.AnnData:
    adata_all = sc.read_10x_h5(h5_path, gex_only=False)
    adata_all.var_names_make_unique()
    is_hto = adata_all.var["feature_types"] == "Multiplexing Capture"
    adata = adata_all[:, ~is_hto].copy()
    annotate_samples(adata)
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
    adata.var["ribo"] = adata.var_names.str.match(r"^(Rpl|Rps)")
    adata.var["hb"] = adata.var_names.str.match(r"^(Hba|Hbb)")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=True, inplace=True)
    return adata


def scenario_mask(
    adata: ad.AnnData,
    count_nmads: float,
    mito_nmads: float,
    max_mito_pct: float | None,
    max_ribo_pct: float | None,
    min_cells_per_sample: int,
) -> tuple[np.ndarray, list[str]]:
    keep = np.ones(adata.n_obs, dtype=bool)
    reasons = []
    for sample_id in adata.obs["sample_id"].unique():
        sample_mask = (adata.obs["sample_id"] == sample_id).to_numpy()
        idx = np.where(sample_mask)[0]
        low_genes, _ = mad_bounds(adata.obs.loc[sample_mask, "log1p_n_genes_by_counts"], count_nmads, "lower")
        low_counts, _ = mad_bounds(adata.obs.loc[sample_mask, "log1p_total_counts"], count_nmads, "lower")
        _, high_mt = mad_bounds(adata.obs.loc[sample_mask, "pct_counts_mt"], mito_nmads, "upper")
        keep[idx] &= adata.obs.loc[sample_mask, "log1p_n_genes_by_counts"].to_numpy() >= low_genes
        keep[idx] &= adata.obs.loc[sample_mask, "log1p_total_counts"].to_numpy() >= low_counts
        keep[idx] &= adata.obs.loc[sample_mask, "pct_counts_mt"].to_numpy() <= high_mt
    if max_mito_pct is not None:
        keep &= adata.obs["pct_counts_mt"].to_numpy() <= max_mito_pct
    if max_ribo_pct is not None:
        keep &= adata.obs["pct_counts_ribo"].to_numpy() <= max_ribo_pct
    kept_by_sample = pd.Series(keep, index=adata.obs["sample_id"]).groupby(level=0).sum()
    low_samples = kept_by_sample[kept_by_sample < min_cells_per_sample]
    if not low_samples.empty:
        reasons.append(f"{len(low_samples)} samples below {min_cells_per_sample} cells")
    return keep, reasons


def plot_distributions(adata: ad.AnnData, out: Path) -> None:
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]
    sc.pl.violin(adata, metrics, groupby="sample_id", rotation=90, show=False)
    plt.savefig(out / "qc_metric_violin_by_sample.png", dpi=150, bbox_inches="tight")
    plt.close()

    for metric in ["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(adata.obs[metric], bins=80)
        ax.set_xlabel(metric)
        ax.set_ylabel("Cells")
        plt.tight_layout()
        plt.savefig(out / f"hist_{metric}.png", dpi=150)
        plt.close()

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(adata.obs["total_counts"], adata.obs["pct_counts_mt"], s=2, alpha=0.25)
    ax.set_xscale("log")
    ax.set_xlabel("total_counts")
    ax.set_ylabel("pct_counts_mt")
    plt.tight_layout()
    plt.savefig(out / "scatter_counts_vs_pct_mt.png", dpi=150)
    plt.close()


def main(h5_path: str, out_dir: str, min_cells_per_sample: int) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = load_qc(h5_path)
    plot_distributions(adata, out)

    by_sample = adata.obs.groupby("sample_id", observed=True).agg(
        n_cells=("sample_id", "size"),
        median_n_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
        p95_pct_mt=("pct_counts_mt", lambda x: np.percentile(x, 95)),
        median_pct_ribo=("pct_counts_ribo", "median"),
        p95_pct_ribo=("pct_counts_ribo", lambda x: np.percentile(x, 95)),
        median_pct_hb=("pct_counts_hb", "median"),
    )
    by_sample.to_csv(out / "qc_metrics_by_sample.csv")

    rows = []
    for count_nmads, mito_nmads, max_mt, max_ribo in product(
        [3.0, 4.0, 5.0, 6.0],
        [2.5, 3.0, 4.0, 5.0],
        [None, 10.0, 15.0, 20.0],
        [None, 60.0, 70.0],
    ):
        keep, reasons = scenario_mask(adata, count_nmads, mito_nmads, max_mt, max_ribo, min_cells_per_sample)
        kept = adata.obs.loc[keep]
        kept_by_sample = kept.groupby("sample_id", observed=True).size()
        rows.append({
            "count_nmads": count_nmads,
            "mito_nmads": mito_nmads,
            "max_mito_pct": "none" if max_mt is None else max_mt,
            "max_ribo_pct": "none" if max_ribo is None else max_ribo,
            "n_cells_keep": int(keep.sum()),
            "n_cells_removed": int((~keep).sum()),
            "pct_removed": float((~keep).mean() * 100),
            "min_cells_per_sample": int(kept_by_sample.min()),
            "median_pct_mt_keep": float(kept["pct_counts_mt"].median()),
            "p95_pct_mt_keep": float(np.percentile(kept["pct_counts_mt"], 95)),
            "median_pct_ribo_keep": float(kept["pct_counts_ribo"].median()),
            "p95_pct_ribo_keep": float(np.percentile(kept["pct_counts_ribo"], 95)),
            "warnings": "; ".join(reasons),
        })
    scan = pd.DataFrame(rows)
    scan = scan.sort_values(["warnings", "pct_removed", "p95_pct_mt_keep"], ascending=[True, True, True])
    scan.to_csv(out / "qc_filter_parameter_scan.csv", index=False)

    print(f"Cells loaded: {adata.n_obs}")
    print(f"Ribosomal genes detected: {int(adata.var['ribo'].sum())}")
    print(f"Hemoglobin genes detected: {int(adata.var['hb'].sum())}")
    print(f"Mitochondrial genes detected: {int(adata.var['mt'].sum())}")
    print(f"Saved QC scan to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-cells-per-sample", type=int, default=100)
    args = parser.parse_args()
    main(args.h5, args.out, args.min_cells_per_sample)
