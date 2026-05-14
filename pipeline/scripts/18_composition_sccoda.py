#!/usr/bin/env python3
"""
18_composition_sccoda.py

Bayesian compositional analysis with scCODA on cell-type frequencies per donor
per condition. Tests:
- Effect of treatment (STM vs DMSO) within each genotype
- Effect of genotype (Mutant vs WT) within each treatment
- Genotype × treatment interaction

Why scCODA: explicit Bayesian model of compositional data, robust to small n,
controls for the simplex constraint (one cell type up → others appear down).
Reference cell type acts as the baseline against which all others are
log-ratio-tested.

Inputs
------
--input      path to adata_progenitor_annotated.h5ad
--out        output directory
--level      'manual_level1' or 'manual_level2'
--ref-type   reference cell type (must be present in all samples; defaults
             to most abundant shared type)

Outputs
-------
sccoda_counts_per_sample.csv         per-donor-per-condition cell-type counts
sccoda_<contrast>_credible_effects.csv  credible effects (FDR-controlled)
sccoda_<contrast>_summary.txt        full HMC sampling summary
sccoda_boxplots_<level>.png          boxplots of celltype fractions × condition
sccoda_chosen_reference.txt          which reference cell type was selected and why
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
    p.add_argument("--out", default="results/18_composition")
    p.add_argument("--level", default="manual_level1",
                   choices=["manual_level1", "manual_level2"])
    p.add_argument("--ref-type", default=None,
                   help="Reference cell type. Defaults to most-abundant type "
                        "with cells in ALL donor×condition samples.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  shape: {adata.shape}, layers: {list(adata.layers.keys())}")

    obs = adata.obs.copy()
    obs["celltype"] = obs[args.level].astype(str)
    obs["condition"] = obs["genotype"].astype(str) + "_" + obs["treatment"].astype(str)
    obs["sample"] = obs["donor"].astype(str) + "__" + obs["condition"].astype(str)

    print("\nDonors per condition:")
    print(obs.groupby("condition")["donor"].nunique())
    print("\nCells per (donor, condition):")
    print(obs.groupby(["condition", "donor"]).size())

    # placeholder for the rest of the workflow
    raise SystemExit("Skeleton OK — Task 3 fills the model")


if __name__ == "__main__":
    main()
