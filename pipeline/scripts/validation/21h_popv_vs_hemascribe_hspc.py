#!/usr/bin/env python3
"""
21h_popv_vs_hemascribe_hspc.py

Within the HSPC compartment (HemaScribe broad == "HSPC"), how does popV's
Nestorowa-vocabulary call (LTHSC / MPP / CMP / GMP / LMPP / MEP) align with
HemaScribe's fine-vocabulary call (HSC / STHSC / MPP2 / MPP3 / MPP4 / MkP /
GMP / CLP)?

Vocab mapping (semantic groups):
  HSC-like:    LTHSC ↔ HSC
  ST-HSC:      (popV has no analogue) ↔ STHSC
  MPP-myelo:   MPP            ↔ MPP3 (FcG+/-), MPP2
  MPP-lymph:   LMPP           ↔ MPP4
  Mk-prog:     MEP/MkP         ↔ MkP
  Myeloid-prog: CMP / GMP     ↔ GMP, NotHSPC
  Lymph-prog:  (CLP)           ↔ CLP

Outputs:
  popv_vs_hemascribe_fine_crosstab.csv      raw counts
  popv_vs_hemascribe_fine_rownorm.csv       row-normalised (popV → HemaScribe distribution)
  popv_vs_hemascribe_summary.md
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default="results/21_hemascribe/validation/popv_vs_hemascribe")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    obs = adata.obs

    in_hspc = obs["hemascribe_broad"] == "HSPC"
    print(f"[21h] HSPC compartment cells: {in_hspc.sum()}")

    sub = obs.loc[in_hspc, ["popv_prediction", "hemascribe_fine",
                             "manual_level2", "population"]].copy()

    # Cross-tab popV vs HemaScribe fine
    ct_counts = pd.crosstab(sub["popv_prediction"], sub["hemascribe_fine"], dropna=False)
    ct_norm = ct_counts.div(ct_counts.sum(axis=1), axis=0).round(3)
    ct_counts.to_csv(os.path.join(args.out, "popv_vs_hemascribe_fine_crosstab.csv"))
    ct_norm.to_csv(os.path.join(args.out, "popv_vs_hemascribe_fine_rownorm.csv"))

    print("\n[21h] Counts:")
    print(ct_counts.to_string())
    print("\n[21h] Row-normalised (popV row → HemaScribe fine distribution):")
    print(ct_norm.to_string())

    # Per-popV-call: top 2 HemaScribe fine matches
    summary_lines = []
    for pop in ct_norm.index:
        s = ct_norm.loc[pop].sort_values(ascending=False)
        s = s[s > 0.01]
        top = ' | '.join(f'{n}({v:.0%})' for n, v in s.head(3).items())
        summary_lines.append(f"  {pop:8s} (n={ct_counts.loc[pop].sum():>5}) -> {top}")

    print("\n[21h] Top HemaScribe fine match per popV call:")
    for ln in summary_lines: print(ln)

    with open(os.path.join(args.out, "popv_vs_hemascribe_summary.md"), "w") as fh:
        fh.write("# popV (Nestorowa) vs HemaScribe fine — HSPC compartment\n\n")
        fh.write(f"Cells: {in_hspc.sum()} (where hemascribe_broad == HSPC)\n\n")
        fh.write("## Top HemaScribe fine matches per popV call\n\n```\n")
        for ln in summary_lines: fh.write(ln + "\n")
        fh.write("\n```\n\n## Row-normalised cross-tab\n\n```\n")
        fh.write(ct_norm.to_string())
        fh.write("\n```\n")

    print(f"\n[21h] Outputs in {args.out}/")


if __name__ == "__main__":
    main()
