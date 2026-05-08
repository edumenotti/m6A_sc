#!/usr/bin/env python
import argparse
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def main(in_path: str, out_dir: str, n_hvgs_values: str) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    values = [int(x) for x in n_hvgs_values.split(",") if x.strip()]
    base = sc.read_h5ad(in_path)
    base.var["mt"] = base.var_names.str.startswith(("mt-", "Mt-", "MT-"))
    base.var["ribo"] = base.var_names.str.match(r"^(Rpl|Rps)")
    base.var["hb"] = base.var_names.str.match(r"^(Hba|Hbb)")

    hvg_sets = {}
    rows = []
    for n_hvgs in values:
        adata = base.copy()
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers["lognorm"] = adata.X.copy()
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvgs,
            flavor="seurat_v3",
            layer="counts",
            batch_key="sample_id",
            span=0.3,
        )
        selected_before = adata.var["highly_variable"].copy()
        n_mt = int((selected_before & adata.var["mt"]).sum())
        n_ribo = int((selected_before & adata.var["ribo"]).sum())
        n_hb = int((selected_before & adata.var["hb"]).sum())
        adata.var.loc[adata.var["mt"] | adata.var["ribo"] | adata.var["hb"], "highly_variable"] = False
        genes = set(adata.var_names[adata.var["highly_variable"]])
        hvg_sets[n_hvgs] = genes

        sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
        var_ratio = adata.uns["pca"]["variance_ratio"]
        rows.append({
            "requested_hvgs": n_hvgs,
            "selected_after_excluding_mt_ribo_hb": len(genes),
            "mt_genes_in_initial_hvgs": n_mt,
            "ribo_genes_in_initial_hvgs": n_ribo,
            "hb_genes_in_initial_hvgs": n_hb,
            "cumvar_pc10": float(var_ratio[:10].sum()),
            "cumvar_pc20": float(var_ratio[:20].sum()),
            "cumvar_pc30": float(var_ratio[:30].sum()),
            "cumvar_pc50": float(var_ratio[:50].sum()),
        })

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(np.arange(1, len(var_ratio) + 1), np.cumsum(var_ratio), marker="o", markersize=2)
        ax.set_xlabel("PC")
        ax.set_ylabel("Cumulative variance ratio")
        ax.set_title(f"{n_hvgs} requested HVGs")
        plt.tight_layout()
        plt.savefig(out / f"pca_cumvar_hvg_{n_hvgs}.png", dpi=150)
        plt.close()

    overlaps = []
    for a, b in combinations(values, 2):
        inter = len(hvg_sets[a] & hvg_sets[b])
        union = len(hvg_sets[a] | hvg_sets[b])
        overlaps.append({
            "hvg_a": a,
            "hvg_b": b,
            "intersection": inter,
            "union": union,
            "jaccard": inter / union if union else np.nan,
        })
    pd.DataFrame(rows).to_csv(out / "hvg_parameter_scan.csv", index=False)
    pd.DataFrame(overlaps).to_csv(out / "hvg_set_overlap.csv", index=False)
    print(f"Saved HVG scan to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--n-hvgs-values", default="1000,2000,3000,4000,5000")
    args = parser.parse_args()
    main(args.input, args.out, args.n_hvgs_values)
