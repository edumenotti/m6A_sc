#!/usr/bin/env python
import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io, sparse


def write_lines(path: Path, values) -> None:
    with path.open("w", newline="") as handle:
        for value in values:
            handle.write(f"{value}\n")


def run_doubletfinder_sample(
    adata,
    sample_id: str,
    helper_r: Path,
    tmp_dir: Path,
    expected_rate: float,
    pcs: int,
    pk: float,
    auto_pk: bool,
    rscript: str,
    install_missing: bool,
    num_cores: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_dir = tmp_dir / sample_id.replace("/", "_")
    sample_dir.mkdir(parents=True, exist_ok=True)
    counts_path = sample_dir / "counts.mtx"
    genes_path = sample_dir / "genes.txt"
    cells_path = sample_dir / "cells.txt"
    calls_path = sample_dir / "doubletfinder_calls.csv"
    summary_path = sample_dir / "doubletfinder_summary.csv"

    counts = adata.X
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    io.mmwrite(counts_path, counts.T.tocoo())
    write_lines(genes_path, adata.var_names)
    write_lines(cells_path, adata.obs_names)

    cmd = [
        rscript,
        str(helper_r),
        "--counts", str(counts_path),
        "--genes", str(genes_path),
        "--cells", str(cells_path),
        "--sample_id", sample_id,
        "--expected_rate", str(expected_rate),
        "--pcs", str(pcs),
        "--pk", str(pk),
        "--num_cores", str(num_cores),
        "--seed", str(seed),
        "--out_csv", str(calls_path),
        "--summary_csv", str(summary_path),
    ]
    if install_missing:
        cmd.append("--install-missing")
    if auto_pk:
        cmd.append("--auto-pk")

    print(f"Running DoubletFinder for {sample_id}: {adata.n_obs} cells")
    subprocess.run(cmd, check=True)
    return pd.read_csv(calls_path), pd.read_csv(summary_path)


def plot_scores(adata, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 4))
    for sample_id in sorted(adata.obs["sample_id"].unique()):
        mask = adata.obs["sample_id"] == sample_id
        scores = adata.obs.loc[mask, "doubletfinder_score"].dropna()
        if scores.empty:
            continue
        ax.hist(scores, bins=50, alpha=0.35, label=sample_id, density=True)
    ax.set_xlabel("DoubletFinder pANN score")
    ax.set_ylabel("Density")
    ax.legend(fontsize=6, ncol=4)
    plt.tight_layout()
    plt.savefig(out / "doublet_scores.png", dpi=150)
    plt.close()


def main(
    in_path: str,
    out_dir: str,
    expected_rate: float,
    pcs: int,
    pk: float,
    auto_pk: bool,
    rscript: str,
    install_missing: bool,
    num_cores: int,
    seed: int,
) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    helper_r = Path(__file__).with_name("run_doubletfinder.R")
    adata = sc.read_h5ad(in_path)

    adata.obs["doubletfinder_score"] = np.nan
    adata.obs["doubletfinder_class"] = "NotRun"
    summaries = []
    calls_all = []

    with tempfile.TemporaryDirectory(prefix="doubletfinder_") as tmp:
        tmp_dir = Path(tmp)
        for sample_id in sorted(adata.obs["sample_id"].unique()):
            mask = adata.obs["sample_id"] == sample_id
            sample = adata[mask].copy()
            calls, summary = run_doubletfinder_sample(
                sample,
                sample_id,
                helper_r,
                tmp_dir,
                expected_rate,
                pcs,
                pk,
                auto_pk,
                rscript,
                install_missing,
                num_cores,
                seed,
            )
            calls["sample_id"] = sample_id
            calls_all.append(calls)
            summaries.append(summary)
            calls_idx = calls.set_index("barcode")
            common = adata.obs_names.intersection(calls_idx.index)
            adata.obs.loc[common, "doubletfinder_score"] = calls_idx.loc[common, "doubletfinder_score"].astype(float)
            adata.obs.loc[common, "doubletfinder_class"] = calls_idx.loc[common, "doubletfinder_class"].astype(str)

    calls_df = pd.concat(calls_all, ignore_index=True)
    summary_df = pd.concat(summaries, ignore_index=True)
    calls_df.to_csv(out / "doubletfinder_calls.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    summary_df.to_csv(out / "doubletfinder_summary.csv", index=False)

    adata.obs["predicted_doublet"] = adata.obs["doubletfinder_class"] == "Doublet"
    plot_scores(adata, out)

    n_doublets = int(adata.obs["predicted_doublet"].sum())
    print(f"Total doublets removed: {n_doublets} ({n_doublets / adata.n_obs * 100:.1f}%)")
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    adata.write_h5ad(out / "adata_no_doublets.h5ad")
    print(f"Saved: {out / 'adata_no_doublets.h5ad'}, {adata.n_obs} cells")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--expected-rate", type=float, default=0.075)
    parser.add_argument("--pcs", type=int, default=20)
    parser.add_argument("--pk", type=float, default=0.09)
    parser.add_argument("--auto-pk", action="store_true")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--install-missing", action="store_true")
    parser.add_argument("--num-cores", type=int, default=1)
    parser.add_argument("--seed", type=int, default=777)
    args = parser.parse_args()
    main(
        args.input,
        args.out,
        args.expected_rate,
        args.pcs,
        args.pk,
        args.auto_pk,
        args.rscript,
        args.install_missing,
        args.num_cores,
        args.seed,
    )
