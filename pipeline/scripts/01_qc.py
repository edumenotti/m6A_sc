#!/usr/bin/env python
import argparse
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

POOL_MAP = {
    "D1_DMSO_LSK": "Pool_A", "D1_DMSO_LK": "Pool_A",
    "D1_STM_LSK": "Pool_B", "D1_STM_LK": "Pool_B",
    "D1_STM_I": "Pool_C", "D1_DMSO_I": "Pool_C",
    "D2_STM_LSK": "Pool_D", "D2_DMSO_LSK": "Pool_D",
    "D2_STM_LK": "Pool_E", "D2_DMSO_LK": "Pool_E",
    "D2_STM_I": "Pool_F", "D2_DMSO_I": "Pool_F",
}


def mad_filter(values: np.ndarray, nmads: float, direction: str = "both") -> np.ndarray:
    """Return a boolean keep mask using MAD thresholds."""
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - med))
    if mad == 0 or np.isnan(mad):
        return np.ones(values.shape, dtype=bool)

    low = med - nmads * mad if direction in ("lower", "both") else -np.inf
    high = med + nmads * mad if direction in ("upper", "both") else np.inf
    return (values >= low) & (values <= high)


def pool_from_sample(sample_id: str) -> str:
    key = "_".join(sample_id.split("_")[:3])
    return POOL_MAP.get(key, "unknown")


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
    meta["pool"] = meta["sample_id"].map(pool_from_sample)
    adata.obs = meta


def plot_qc(adata: ad.AnnData, out_path: Path) -> None:
    keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    if "pct_counts_ribo" in adata.obs:
        keys.append("pct_counts_ribo")
    if "pct_counts_hb" in adata.obs:
        keys.append("pct_counts_hb")
    sc.pl.violin(
        adata,
        keys,
        groupby="sample_id",
        rotation=90,
        show=False,
    )
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def main(h5_path: str, out_dir: str, mito_nmads: float, count_nmads: float) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    adata_all = sc.read_10x_h5(h5_path, gex_only=False)
    adata_all.var_names_make_unique()

    is_hto = adata_all.var["feature_types"] == "Multiplexing Capture"
    hto = adata_all[:, is_hto].copy()
    adata = adata_all[:, ~is_hto].copy()
    adata.obsm["HTO"] = hto.X.toarray() if hasattr(hto.X, "toarray") else np.asarray(hto.X)
    adata.uns["HTO_names"] = hto.var_names.tolist()

    annotate_samples(adata)

    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
    adata.var["ribo"] = adata.var_names.str.match(r"^(Rpl|Rps)")
    adata.var["hb"] = adata.var_names.str.match(r"^(Hba|Hbb)")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=True, inplace=True)

    keep = np.ones(adata.n_obs, dtype=bool)
    for sample_id in adata.obs["sample_id"].unique():
        mask = (adata.obs["sample_id"] == sample_id).to_numpy()
        idx = np.where(mask)[0]
        keep[idx] &= mad_filter(adata.obs.loc[mask, "log1p_n_genes_by_counts"], count_nmads, "lower")
        keep[idx] &= mad_filter(adata.obs.loc[mask, "log1p_total_counts"], count_nmads, "lower")
        keep[idx] &= mad_filter(adata.obs.loc[mask, "pct_counts_mt"], mito_nmads, "upper")

    before_summary = adata.obs.groupby("sample_id", observed=True).size().rename("n_cells_before")
    plot_qc(adata, out / "qc_violin_before.png")

    print(f"Cells before QC: {adata.n_obs}")
    print(f"Cells removed: {int((~keep).sum())}")
    adata = adata[keep].copy()
    print(f"Cells after QC: {adata.n_obs}")

    plot_qc(adata, out / "qc_violin_after.png")

    after_summary = adata.obs.groupby("sample_id", observed=True).size().rename("n_cells_after")
    summary = pd.concat([before_summary, after_summary], axis=1).fillna(0).astype(int).reset_index()
    summary["n_removed"] = summary["n_cells_before"] - summary["n_cells_after"]
    metric_summary = adata.obs.groupby("sample_id", observed=True).agg(
        median_n_genes=("n_genes_by_counts", "median"),
        median_total_counts=("total_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
        median_pct_ribo=("pct_counts_ribo", "median"),
        median_pct_hb=("pct_counts_hb", "median"),
    ).reset_index()
    summary = summary.merge(metric_summary, on="sample_id", how="left")
    summary.to_csv(out / "cells_per_sample.csv", index=False)

    adata.write_h5ad(out / "adata_qc.h5ad")
    print(f"Saved: {out / 'adata_qc.h5ad'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--mito_nmads", type=float, default=4.0)
    parser.add_argument("--count_nmads", type=float, default=6.0)
    args = parser.parse_args()
    main(args.h5, args.out, args.mito_nmads, args.count_nmads)
