#!/usr/bin/env python3
"""
21a_export_for_hemascribe.py

Export AnnData counts + obs metadata to a 10x-compatible MTX bundle that
Seurat can read directly via Seurat::ReadMtx. This sidesteps zellkonverter
(which tries to bootstrap its own Python via basilisk/reticulate and fails
under pixi).

Outputs (in --out dir):
  counts.mtx.gz
  barcodes.tsv.gz
  features.tsv.gz
  obs.csv               full obs metadata (used by R to attach @meta.data)

Usage:
  pixi run python pipeline/scripts/21a_export_for_hemascribe.py \\
    --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \\
    --out   results/21_hemascribe/export \\
    --layer counts
"""
import argparse
import gzip
import os
import shutil

import scanpy as sc
import scipy.io
import scipy.sparse


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--layer", default="counts")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"[21a] Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"      shape: {adata.shape}")

    if args.layer not in adata.layers:
        raise SystemExit(f"layer '{args.layer}' not in adata.layers={list(adata.layers)}")

    X = adata.layers[args.layer]
    if not scipy.sparse.isspmatrix(X):
        X = scipy.sparse.csr_matrix(X)
    # MTX is gene-x-cell for the 10x convention; AnnData is cell-x-gene.
    X_gxc = X.T.tocoo()

    mtx_path = os.path.join(args.out, "counts.mtx")
    print(f"[21a] Writing {mtx_path} ({X_gxc.shape[0]} features x {X_gxc.shape[1]} cells)")
    scipy.io.mmwrite(mtx_path, X_gxc, field="integer")
    with open(mtx_path, "rb") as fin, gzip.open(mtx_path + ".gz", "wb") as fout:
        shutil.copyfileobj(fin, fout)
    os.remove(mtx_path)

    bc_path = os.path.join(args.out, "barcodes.tsv.gz")
    with gzip.open(bc_path, "wt") as fh:
        for b in adata.obs_names:
            fh.write(f"{b}\n")
    print(f"[21a] Wrote {bc_path}")

    # features.tsv: ID, name, type — Seurat::ReadMtx feature.column=1 uses col1
    ft_path = os.path.join(args.out, "features.tsv.gz")
    with gzip.open(ft_path, "wt") as fh:
        for g in adata.var_names:
            fh.write(f"{g}\t{g}\tGene Expression\n")
    print(f"[21a] Wrote {ft_path}")

    obs_path = os.path.join(args.out, "obs.csv")
    adata.obs.to_csv(obs_path)
    print(f"[21a] Wrote {obs_path}")

    print("[21a] Export complete.")


if __name__ == "__main__":
    main()
