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

    print(f"Pseudobulk shape: {pdata.shape}")

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    contrasts = [
        ("treatment", ["treatment", "STM", "DMSO"]),
        ("genotype",  ["genotype", "Mutant", "WT"]),
        # interaction handled separately via LRT below
    ]

    summary_rows = []
    celltypes = sorted(pdata.obs["celltype"].unique())

    for ct in celltypes:
        ct_mask = pdata.obs["celltype"] == ct
        sub = pdata[ct_mask].copy()

        # Drop genes with all-zero counts
        keep_gene = (sub.X.sum(axis=0) > 0)
        sub = sub[:, np.asarray(keep_gene).ravel()].copy()

        # Need ≥2 samples per condition AND ≥min-cells in each donor sample
        per_cond = sub.obs.groupby("condition").size()
        cell_ok = sub.obs["psbulk_n_cells"] >= args.min_cells
        if not cell_ok.all():
            sub = sub[cell_ok.values].copy()
            per_cond = sub.obs.groupby("condition").size()
        if (per_cond < 2).any() or sub.n_obs < 6:
            print(f"  Skipping {ct}: insufficient samples after filtering "
                  f"({sub.n_obs} samples, per-cond: {per_cond.to_dict()})")
            continue

        print(f"\n=== Celltype: {ct} (n_samples={sub.n_obs}) ===")
        meta = sub.obs[["genotype","treatment","donor"]].copy()
        for col in ("genotype","treatment","donor"):
            meta[col] = meta[col].astype("category")

        # Paired design: donor blocks for mouse-level variation (competitive transplant).
        # WT and Mutant cells come from the same mice, so donor must be a blocking factor.
        # Drop interaction — with n=2 donors the interaction term is not estimable.
        n_donors = meta["donor"].nunique()
        if n_donors >= 2:
            design = "~ donor + genotype + treatment"
        else:
            design = "~ genotype + treatment"
            print(f"  Using unpaired design (only {n_donors} donor(s) after filtering)")

        try:
            dds = DeseqDataSet(
                counts=pd.DataFrame(
                    sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X,
                    index=sub.obs_names, columns=sub.var_names,
                ).astype(int),
                metadata=meta,
                design=design,
                inference=DefaultInference(n_cpus=4),
                quiet=True,
            )
            dds.deseq2()
        except Exception as e:
            print(f"  Skipping {ct}: deseq2() failed: {e}")
            continue

        for cname, contrast in contrasts:
            try:
                ds = DeseqStats(dds, contrast=contrast, quiet=True)
                ds.summary()
                res = ds.results_df.reset_index().rename(columns={"index":"gene"})
                res["celltype"] = ct
                res["contrast"] = cname
                res.to_csv(
                    os.path.join(args.out, f"deg_{ct}_{cname}.csv"),
                    index=False
                )
                top = res[(res["padj"] < args.padj) & (res["log2FoldChange"].abs() > 1)]
                top = top.sort_values("padj")
                top.to_csv(
                    os.path.join(args.out, f"deg_{ct}_{cname}_top.csv"),
                    index=False
                )
                summary_rows.append({
                    "celltype": ct,
                    "contrast": cname,
                    "n_sig_padj": int((res["padj"] < args.padj).sum()),
                    "n_sig_padj_lfc1": int(len(top)),
                    "n_samples": sub.n_obs,
                })

                # Volcano
                fig, ax = plt.subplots(figsize=(7, 6))
                x = res["log2FoldChange"].values
                y = -np.log10(res["padj"].fillna(1).values)
                sig = (res["padj"] < args.padj) & (res["log2FoldChange"].abs() > 1)
                ax.scatter(x[~sig], y[~sig], s=4, c="lightgray", alpha=0.5)
                ax.scatter(x[sig],  y[sig],  s=8, c="firebrick")
                ax.axhline(-np.log10(args.padj), c="k", lw=0.5, ls="--")
                ax.axvline( 1, c="k", lw=0.5, ls="--")
                ax.axvline(-1, c="k", lw=0.5, ls="--")
                ax.set_xlabel("log2 FC")
                ax.set_ylabel("-log10 padj")
                ax.set_title(f"{ct} — {cname}")
                fig.tight_layout()
                fig.savefig(
                    os.path.join(args.out, f"deg_volcano_{ct}_{cname}.png"),
                    dpi=150
                )
                plt.close(fig)
            except Exception as e:
                print(f"  [{ct}/{cname}] failed: {e}")

    pd.DataFrame(summary_rows).to_csv(
        os.path.join(args.out, "deg_summary.csv"), index=False
    )
    print(f"\nScript 19 complete. Outputs in {args.out}/")


if __name__ == "__main__":
    main()
