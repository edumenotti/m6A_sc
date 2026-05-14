#!/usr/bin/env python3
"""
19_pseudobulk_deg.py

Per-celltype pseudobulk differential expression with PyDESeq2.
Donor is the replicate unit (n=2 per condition → 8 donor-level samples per celltype).

For each celltype with ≥10 cells in every (donor × condition) sample, we:
  1. Sum raw counts across cells → (donor × condition) × gene matrix
  2. Fit DESeq2 with design ~genotype + treatment + genotype:treatment
  3. Extract Wald results for treatment, genotype, and interaction
  4. Write per-celltype CSVs + a global summary

Why this is appropriate:
- Donor as replicate avoids pseudoreplication (the CellChat sin).
- DESeq2 dispersion shrinkage stabilises n=2 variance estimates.
- Restricting to celltypes present in every sample removes the composition
  confound that broke CellChat.

Inputs
------
--input        path to adata_progenitor_annotated.h5ad
--out          output directory
--level        manual_level1 (default) or manual_level2
--min-cells    min cells per (donor, condition, celltype) to retain (default 10)
--padj         FDR threshold for "significant" gene tables (default 0.05)

Outputs
-------
pseudobulk_samples_per_celltype.csv     QC: cells/sample/celltype matrix
deg_<celltype>_<contrast>.csv           full DESeq2 results
deg_<celltype>_<contrast>_top.csv       top sig hits (|lfc|>1, padj<thresh)
deg_summary.csv                         # of sig genes per celltype × contrast
deg_volcano_<celltype>_<contrast>.png   volcano plots per (celltype, contrast)
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default="results/19_pseudobulk_deg")
    p.add_argument("--level", default="manual_level1",
                   choices=["manual_level1", "manual_level2"])
    p.add_argument("--min-cells", type=int, default=10)
    p.add_argument("--padj", type=float, default=0.05)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    assert "counts" in adata.layers, "Need adata.layers['counts'] (raw integer counts)"

    adata.obs["celltype"] = adata.obs[args.level].astype(str)
    adata.obs["condition"] = (
        adata.obs["genotype"].astype(str) + "_" + adata.obs["treatment"].astype(str)
    )
    adata.obs["sample"] = (
        adata.obs["donor"].astype(str) + "__" + adata.obs["condition"].astype(str)
    )

    # ── Pseudobulk via decoupler (sum, raw counts layer) ────────────────
    import decoupler as dc
    pdata = dc.pp.pseudobulk(
        adata,
        sample_col="sample",
        groups_col="celltype",
        layer="counts",
        mode="sum",
        skip_checks=False,
    )
    pdata.obs["donor"]     = pdata.obs["sample"].str.split("__").str[0]
    pdata.obs["condition"] = pdata.obs["sample"].str.split("__").str[1]
    pdata.obs["genotype"]  = pdata.obs["condition"].str.split("_").str[0]
    pdata.obs["treatment"] = pdata.obs["condition"].str.split("_").str[1]

    # QC table: cells per (sample, celltype)
    # decoupler v2 renamed the cell-count column from 'psbulk_n_cells' to 'psbulk_cells'.
    n_cells_col = "psbulk_cells" if "psbulk_cells" in pdata.obs.columns else "psbulk_n_cells"
    pdata.obs["psbulk_n_cells"] = pdata.obs[n_cells_col]
    qc = pdata.obs[["sample","celltype","psbulk_n_cells"]].copy()
    qc.to_csv(os.path.join(args.out, "pseudobulk_samples_per_celltype.csv"), index=False)

    # placeholder for DEG loop (Task 5)
    print(f"Pseudobulk shape: {pdata.shape}")
    raise SystemExit("Aggregation OK — Task 5 adds PyDESeq2")


if __name__ == "__main__":
    main()
