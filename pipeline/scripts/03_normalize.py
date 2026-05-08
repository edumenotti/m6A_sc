#!/usr/bin/env python
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


def main(in_path: str, out_dir: str, n_hvgs: int) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
    adata.var["ribo"] = adata.var_names.str.match(r"^(Rpl|Rps)")
    adata.var["hb"] = adata.var_names.str.match(r"^(Hba|Hbb)")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvgs,
            flavor="seurat_v3",
            layer="counts",
            batch_key="sample_id",
            span=0.3,
        )
    except ImportError:
        print("scikit-misc is unavailable; falling back to cell_ranger HVG flavor.")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvgs,
            flavor="cell_ranger",
            batch_key="sample_id",
        )

    selected_before = adata.var["highly_variable"].copy()
    excluded = adata.var["mt"] | adata.var["ribo"] | adata.var["hb"]
    adata.var.loc[excluded, "highly_variable"] = False
    print(
        "HVGs selected: "
        f"{int(adata.var['highly_variable'].sum())} "
        f"(excluded mt={int((selected_before & adata.var['mt']).sum())}, "
        f"ribo={int((selected_before & adata.var['ribo']).sum())}, "
        f"hb={int((selected_before & adata.var['hb']).sum())})"
    )

    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(out / "hvg_dispersion.png", dpi=150, bbox_inches="tight")
    plt.close()

    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False)
    plt.savefig(out / "pca_variance_ratio.png", dpi=150, bbox_inches="tight")
    plt.close()

    adata.write_h5ad(out / "adata_normalized.h5ad")
    print(f"Saved: {out / 'adata_normalized.h5ad'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--n_hvgs", type=int, default=3000)
    args = parser.parse_args()
    main(args.input, args.out, args.n_hvgs)
