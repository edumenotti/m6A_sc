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

    # ── Build sample × celltype count matrix ─────────────────────────────
    counts = (
        obs.groupby(["sample", "celltype"]).size().unstack(fill_value=0).reset_index()
    )
    # Add covariates from the sample name (donor + condition)
    counts["donor"] = counts["sample"].str.split("__").str[0]
    counts["condition"] = counts["sample"].str.split("__").str[1]
    counts["genotype"] = counts["condition"].str.split("_").str[0]
    counts["treatment"] = counts["condition"].str.split("_").str[1]
    counts.to_csv(os.path.join(args.out, "sccoda_counts_per_sample.csv"), index=False)
    print(f"\nCount matrix shape: {counts.shape}")
    print(counts.head().to_string())

    # ── Reference cell type selection ────────────────────────────────────
    celltype_cols = [c for c in counts.columns
                     if c not in ("sample", "donor", "condition", "genotype", "treatment")]
    # Most-abundant cell type with >0 in every sample
    presence = (counts[celltype_cols] > 0).all(axis=0)
    candidates = [c for c in celltype_cols if presence[c]]
    if args.ref_type:
        ref = args.ref_type
        assert ref in celltype_cols, f"--ref-type {ref} not in {celltype_cols}"
    elif candidates:
        totals = counts[candidates].sum(axis=0).sort_values(ascending=False)
        ref = totals.index[0]
    else:
        # fallback: use scCODA's automatic selector via 'automatic'
        ref = "automatic"
    with open(os.path.join(args.out, "sccoda_chosen_reference.txt"), "w") as fh:
        fh.write(f"reference_cell_type: {ref}\n")
        fh.write(f"present_in_all_samples: {candidates}\n")
        fh.write(f"sample_totals_per_celltype:\n{counts[celltype_cols].sum(axis=0).to_string()}\n")
    print(f"\nReference cell type: {ref}")

    # ── scCODA model: factorial formula ──────────────────────────────────
    from sccoda.util import cell_composition_data as dat
    from sccoda.util import comp_ana as mod

    sccoda_data = dat.from_pandas(counts, covariate_columns=[
        "sample", "donor", "condition", "genotype", "treatment"
    ])

    contrasts = [
        ("treatment_in_WT",         "treatment", counts["genotype"] == "WT"),
        ("treatment_in_Mutant",     "treatment", counts["genotype"] == "Mutant"),
        ("genotype_in_DMSO",        "genotype",  counts["treatment"] == "DMSO"),
        ("genotype_in_STM",         "genotype",  counts["treatment"] == "STM"),
        ("genotype_x_treatment",    "genotype + treatment + genotype:treatment", None),
    ]

    for name, formula, row_mask in contrasts:
        print(f"\n=== Contrast: {name} (formula: {formula}) ===")
        if row_mask is not None:
            sub = sccoda_data[row_mask.values].copy()
        else:
            sub = sccoda_data.copy()
        if sub.shape[0] < 2:
            print(f"  skipping — only {sub.shape[0]} samples after filtering")
            continue
        model = mod.CompositionalAnalysis(
            sub, formula=formula, reference_cell_type=ref
        )
        result = model.sample_hmc(num_results=20000, num_burnin=5000)
        result.set_fdr(est_fdr=0.05)
        ce = result.credible_effects()
        ce_df = ce.reset_index()
        ce_df.columns = ["covariate", "celltype", "credible"]
        ce_df.to_csv(
            os.path.join(args.out, f"sccoda_{name}_credible_effects.csv"),
            index=False
        )
        with open(os.path.join(args.out, f"sccoda_{name}_summary.txt"), "w") as fh:
            import io, contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                result.summary_extended()
            fh.write(buf.getvalue())
        print(f"  credible effects:\n{ce_df[ce_df['credible']].to_string(index=False)}")

    # ── Diagnostic boxplots ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    long = counts.melt(
        id_vars=["sample", "donor", "condition", "genotype", "treatment"],
        value_vars=celltype_cols,
        var_name="celltype", value_name="n_cells",
    )
    long["fraction"] = long.groupby("sample")["n_cells"].transform(
        lambda x: x / x.sum()
    )
    long.boxplot(column="fraction", by=["celltype", "condition"], ax=axes[0])
    axes[0].set_xticklabels([])
    axes[0].set_title(f"Fractions per condition ({args.level})")
    axes[0].set_xlabel("celltype × condition")
    plt.suptitle("")
    pivot = (long.groupby(["celltype", "condition"])["fraction"].mean()
             .unstack(fill_value=0))
    pivot.plot(kind="bar", stacked=True, ax=axes[1])
    axes[1].set_title("Mean fraction per condition")
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fig.savefig(os.path.join(args.out, f"sccoda_boxplots_{args.level}.png"), dpi=150)
    plt.close(fig)

    print(f"\nscCODA analysis complete. Outputs in {args.out}/")


if __name__ == "__main__":
    main()
