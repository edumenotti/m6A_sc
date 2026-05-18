#!/usr/bin/env python3
"""
21g_stem_gate_outliers.py

Triage the 303 cells in the LSK (Lin-/Sca1+/c-Kit+) sort gate that HemaScribe
calls mature lineage or committed progenitor — are they (a) real Sca1+ cells
caught in a transition (HemaScribe sees committed transcriptome but Sca1 RNA
is still on), or (b) sort contamination (Sca1 negative on RNA → leak through
gate)?

For each suspect subgroup, report mean Sca1 (Ly6a) and Procr (HSC) expression.
Compare against the LSK "core" cells (HemaScribe HSPC + manual HSC/MPP).

Subgroups inspected:
  - LSK + HemaScribe Plasma_cell      (31 — strong leak suspect)
  - LSK + HemaScribe cMoP             (87)
  - LSK + HemaScribe mGMP             (58)
  - LSK + HemaScribe EryP             (45)
  - LSK + HemaScribe other-mature     (~50)
  - reference: LSK + HemaScribe HSPC  (the canonical bulk)

Outputs:
  lsk_outliers_marker_table.csv
  lsk_outliers_dotplot.png
  lsk_outliers_verdict.md
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def to_dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)


MARKERS = {
    "Stem_gate":  ["Ly6a", "Kit", "Cd34"],   # gate-defining markers
    "HSC":        ["Procr", "Hlf", "Mecom"],
    "Erythroid":  ["Klf1", "Gata1", "Hba-a1"],
    "Myeloid":    ["Mpo", "Elane", "Prtn3", "Cebpa"],
    "Plasma":     ["Jchain", "Mzb1", "Xbp1", "Sdc1"],
    "Mature_B":   ["Ighm", "Cd19", "Ms4a1"],
}
ALL = sorted({g for v in MARKERS.values() for g in v})


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default="results/21_hemascribe/validation/lsk_outliers")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"]

    available = [g for g in ALL if g in adata.var_names]
    missing = sorted(set(ALL) - set(available))
    if missing:
        print(f"[21g] missing genes: {missing}")

    obs = adata.obs.copy()
    sub = pd.Series("other", index=obs.index, dtype=object)

    in_lsk = obs["population"] == "LSK"
    sub[in_lsk & (obs["hemascribe_broad"] == "Plasma_cell")]   = "LSK_Plasma"
    sub[in_lsk & (obs["hemascribe_broad"] == "cMoP")]          = "LSK_cMoP"
    sub[in_lsk & (obs["hemascribe_broad"] == "mGMP")]          = "LSK_mGMP"
    sub[in_lsk & (obs["hemascribe_broad"] == "EryP")]          = "LSK_EryP"
    sub[in_lsk & obs["hemascribe_broad"].isin(
        ["Neutrophil", "Monocyte", "B_cell", "Immature_B_cell",
         "T_cell", "Mast_cell", "RBC", "cDC", "Megakaryocyte", "GP"])] = "LSK_other_mature"
    sub[in_lsk & (obs["hemascribe_broad"] == "HSPC")]          = "LSK_HSPC_reference"
    adata.obs["lsk_subgroup"] = pd.Categorical(sub)

    targets = ["LSK_Plasma", "LSK_cMoP", "LSK_mGMP", "LSK_EryP",
               "LSK_other_mature", "LSK_HSPC_reference"]
    print("\n[21g] subgroup sizes:")
    print(adata.obs["lsk_subgroup"].value_counts().reindex(targets, fill_value=0).to_string())

    rows = []
    for grp in targets:
        cells = adata.obs.index[adata.obs["lsk_subgroup"] == grp]
        if len(cells) == 0: continue
        row = {"subgroup": grp, "n_cells": len(cells)}
        for g in available:
            x = to_dense(adata[cells, g].X).ravel()
            row[f"{g}_mean"]    = float(x.mean())
            row[f"{g}_pct_pos"] = float((x > 0).mean() * 100)
        rows.append(row)
    table = pd.DataFrame(rows).set_index("subgroup")
    table.to_csv(os.path.join(args.out, "lsk_outliers_marker_table.csv"))

    pct = [c for c in table.columns if c.endswith("_pct_pos")]
    print("\n[21g] %-positive per marker:")
    print(table[pct].round(0).to_string())

    sc.pl.dotplot(adata, var_names={k: [g for g in v if g in available]
                                     for k, v in MARKERS.items()},
                  groupby="lsk_subgroup", standard_scale="var", show=False)
    plt.savefig(os.path.join(args.out, "lsk_outliers_dotplot.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    with open(os.path.join(args.out, "lsk_outliers_verdict.md"), "w") as fh:
        fh.write("# LSK gate outlier triage\n\n")
        fh.write("If a subgroup has Sca1 (Ly6a) RNA mean comparable to LSK_HSPC_reference,\n")
        fh.write("the cell really was Sca1+ at sort time (transition / Sca1 RNA persistence).\n")
        fh.write("If Sca1 mean is much lower than reference, sort contamination is the\n")
        fh.write("more likely explanation.\n\n")
        fh.write("Plasma_cell in LSK is almost certainly leak (mature B-derived plasma\n")
        fh.write("cells should be Lin+).\n\n")
        fh.write("## %-positive per marker per subgroup\n\n```\n")
        fh.write(table[pct].round(0).to_string())
        fh.write("\n```\n")
    print(f"[21g] Outputs in {args.out}/")


if __name__ == "__main__":
    main()
