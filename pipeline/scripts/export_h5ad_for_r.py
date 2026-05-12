#!/usr/bin/env python3
"""
export_h5ad_for_r.py

Export AnnData → MatrixMarket bundle for R-side analyses (CellChat, NicheNet).
Avoids zellkonverter/basilisk's pyenv-build path on Ubuntu 24.04.

Exports RAW COUNTS by default (from adata.layers["counts"]) since downstream
R packages (CellChat, NicheNet) call their own normalizeData() and assume
input is unnormalized. Pass --layer lognorm or --layer X to override.

Inputs
------
--input  : path to .h5ad
--out    : output directory
--layer  : which expression matrix to export. One of:
             "counts"  (default — adata.layers["counts"], raw integer)
             "lognorm" (adata.layers["lognorm"], log1p-normalized)
             "X"       (adata.X — whatever is currently there)

Outputs (in --out/)
-------
counts.mtx     : sparse matrix, genes × cells (CellChat convention)
barcodes.txt   : cell barcodes (one per line)
features.csv   : gene metadata (var)
metadata.csv   : cell metadata (obs)
export_info.json : layer name, integer-ness, shape (for sanity-check from R)
"""
import argparse
import json
import os
import warnings

import numpy as np
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

warnings.filterwarnings("ignore", category=FutureWarning)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default=".")
    p.add_argument("--layer", default="counts", choices=["counts", "lognorm", "X"])
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  shape: {adata.shape}, layers: {list(adata.layers.keys())}")

    if args.layer == "X":
        mat = adata.X
    elif args.layer not in adata.layers:
        raise SystemExit(f"layer '{args.layer}' not in adata.layers: {list(adata.layers.keys())}")
    else:
        mat = adata.layers[args.layer]

    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)

    # CellChat expects genes-as-rows
    mat_gxc = mat.T.tocsr()

    sample = mat_gxc.data[: min(10000, mat_gxc.data.size)]
    is_integer = bool(np.all(sample == sample.astype(np.int64)))
    print(f"  exporting layer='{args.layer}' shape={mat_gxc.shape} (genes×cells), integer={is_integer}")

    mtx_path = os.path.join(args.out, "counts.mtx")
    print(f"  writing {mtx_path} ...")
    sio.mmwrite(mtx_path, mat_gxc, field="integer" if is_integer else "real")

    with open(os.path.join(args.out, "barcodes.txt"), "w") as fh:
        fh.write("\n".join(adata.obs_names.astype(str)) + "\n")

    adata.var.to_csv(os.path.join(args.out, "features.csv"))
    adata.obs.to_csv(os.path.join(args.out, "metadata.csv"))

    with open(os.path.join(args.out, "export_info.json"), "w") as fh:
        json.dump(
            {
                "layer": args.layer,
                "n_genes": int(mat_gxc.shape[0]),
                "n_cells": int(mat_gxc.shape[1]),
                "integer": is_integer,
                "orientation": "genes_rows_cells_cols",
            },
            fh,
            indent=2,
        )

    print(f"Done. Outputs in {args.out}/")


if __name__ == "__main__":
    main()
