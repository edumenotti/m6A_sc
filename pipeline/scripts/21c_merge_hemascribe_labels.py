#!/usr/bin/env python3
"""
21c_merge_hemascribe_labels.py

Companion to 21b_hemascribe.R. Reads the HemaScribe label CSV, merges into
the AnnData object, applies sort-gate biological priors (LSK / LK / I) to
flag/reject implausible calls, and writes audit outputs:

  - results/14_progenitor_annotated/adata_hemascribe.h5ad
      with new obs columns:
        hemascribe_broad
        hemascribe_fine
        hemascribe_score
        hemascribe_broad_sortgated      (NaN where sort-gate rejects)
        hemascribe_fine_sortgated       (NaN where sort-gate rejects)
        hemascribe_sortgate_flag        (True where any call violates prior)

  - results/21_hemascribe/
        confusion_manual_level1_vs_hemascribe_broad.csv
        confusion_manual_level2_vs_hemascribe_fine.csv
        confusion_sort_vs_hemascribe_broad_raw.csv
        confusion_sort_vs_hemascribe_fine_raw.csv
        sortgate_flag_rates.csv
        confusion_heatmap_level2.png
        README.md

Usage:
  pixi run python pipeline/scripts/21c_merge_hemascribe_labels.py \\
    --adata results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \\
    --labels results/21_hemascribe/hemascribe_labels.csv \\
    --out results/14_progenitor_annotated/adata_hemascribe.h5ad \\
    --audit-dir results/21_hemascribe
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


# Sort-gate priors. Labels are LOWERCASED before comparison so we tolerate
# variations like "HSC" / "hsc" / "Hsc". The check is substring-permissive:
# any allowed token appearing in the label passes (covers "LT-HSC", "MPP4",
# "GMP_neutrophil_primed" etc.).
SORT_PRIORS = {
    # LSK = Lin-/Sca1+/c-Kit+ → stem / multipotent / early lymphoid (CLP is Sca1+).
    # HemaScribe labels seen here: HSPC (broad), HSC/STHSC/MPP{2,3,4}/FcG_*_MPP3 (fine).
    "LSK": {
        "hsc", "mpp", "lmpp", "lsk", "stemcell", "stem", "clp",
        "hspc", "sthsc",
    },
    # LK  = Lin-/Sca1-/c-Kit+ → committed myeloid/erythroid progenitors.
    # HemaScribe broad here: GP, mGMP, EryP, cMoP, MkP, Megakaryocyte (and the fine GMP/MkP/CLP).
    "LK": {
        "cmp", "gmp", "mep", "mkp", "emp", "cfu", "myeloid", "erythroid",
        "megakaryocyte", "progenitor", "granulocyte",
        "eryp", "cmop", "gp", "mgmp",
    },
    # I   = ungated / broad → all labels allowed
    "I": None,
}


# Labels HemaScribe uses as explicit "no call" sentinels (not a misclassification —
# the classifier simply did not assign that granularity to the cell, e.g. NotHSPC
# is the default for cells outside the HSPC compartment when fine_classify runs).
HEMA_NO_CALL = {"nothspc", "nothem", "nan", "na", ""}


def label_passes_prior(label: str, allowed_tokens):
    if allowed_tokens is None:
        return True
    if not isinstance(label, str) or not label:
        return True  # NaN/empty: don't flag
    s = label.lower()
    if s in HEMA_NO_CALL:
        return True  # deliberate no-call by HemaScribe; treat like NaN
    return any(tok in s for tok in allowed_tokens)


def apply_sort_priors(labels: pd.Series, sort: pd.Series) -> pd.DataFrame:
    """Return df with raw label, sortgated label (NaN where rejected), and flag."""
    raw = labels.copy()
    gated = labels.astype(object).copy()
    flag = pd.Series(False, index=labels.index)
    for gate, allowed in SORT_PRIORS.items():
        mask = (sort == gate)
        if not mask.any() or allowed is None:
            continue
        bad = mask & ~labels.fillna("").map(lambda x: label_passes_prior(x, allowed))
        flag[bad] = True
        gated[bad] = np.nan
    return pd.DataFrame({
        "raw": raw,
        "sortgated": gated,
        "flag": flag,
    })


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True)
    p.add_argument("--labels", required=True,
                   help="CSV from 21b_hemascribe.R with cell_id, hemascribe_broad, hemascribe_fine, hemascribe_score")
    p.add_argument("--out", required=True)
    p.add_argument("--audit-dir", default="results/21_hemascribe")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.audit_dir, exist_ok=True)
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f"[21c] Loading {args.adata}")
    adata = sc.read_h5ad(args.adata)
    print(f"      shape: {adata.shape}")

    print(f"[21c] Loading labels {args.labels}")
    lab = pd.read_csv(args.labels)
    lab = lab.set_index("cell_id")
    lab = lab.reindex(adata.obs_names)  # AnnData order; missing cells get NaN
    n_missing = (lab["hemascribe_broad"].isna() | lab["hemascribe_fine"].isna()).sum()
    print(f"      cells with NaN in broad or fine: {n_missing}")

    adata.obs["hemascribe_broad"] = pd.Categorical(lab["hemascribe_broad"])
    adata.obs["hemascribe_fine"]  = pd.Categorical(lab["hemascribe_fine"])
    adata.obs["hemascribe_score"] = lab["hemascribe_score"].values

    # ── Sort-gate priors on both granularities ──────────────────────────
    sort = adata.obs["population"].astype(str)

    broad_audit = apply_sort_priors(lab["hemascribe_broad"], sort)
    fine_audit  = apply_sort_priors(lab["hemascribe_fine"],  sort)

    adata.obs["hemascribe_broad_sortgated"] = pd.Categorical(broad_audit["sortgated"])
    adata.obs["hemascribe_fine_sortgated"]  = pd.Categorical(fine_audit["sortgated"])
    adata.obs["hemascribe_sortgate_flag"]   = (broad_audit["flag"] | fine_audit["flag"]).values

    n = adata.n_obs
    n_flag_broad = int(broad_audit["flag"].sum())
    n_flag_fine  = int(fine_audit["flag"].sum())
    n_flag_any   = int((broad_audit["flag"] | fine_audit["flag"]).sum())
    print(f"\n[21c] Sort-gate conflicts:")
    print(f"      broad: {n_flag_broad}/{n} ({n_flag_broad/n*100:.1f}%)")
    print(f"      fine : {n_flag_fine}/{n} ({n_flag_fine/n*100:.1f}%)")
    print(f"      any  : {n_flag_any}/{n} ({n_flag_any/n*100:.1f}%)")

    flag_rates = adata.obs.groupby("population", observed=True).apply(
        lambda d: pd.Series({
            "n_cells": len(d),
            "flag_broad_pct": 100 * broad_audit.loc[d.index, "flag"].mean(),
            "flag_fine_pct":  100 * fine_audit.loc[d.index, "flag"].mean(),
        })
    )
    print("\nFlag rate per sort gate:")
    print(flag_rates.to_string())
    flag_rates.to_csv(os.path.join(args.audit_dir, "sortgate_flag_rates.csv"))

    # ── Confusion matrices ──────────────────────────────────────────────
    def crosstab_save(rows_col, cols_col, name, normalize="index"):
        ct = pd.crosstab(
            adata.obs[rows_col].astype(str),
            adata.obs[cols_col].astype(str),
            dropna=False,
            normalize=normalize if normalize else False,
        )
        ct.to_csv(os.path.join(args.audit_dir, f"{name}.csv"))
        return ct

    crosstab_save("manual_level1", "hemascribe_broad_sortgated",
                  "confusion_manual_level1_vs_hemascribe_broad")
    conf_fine = crosstab_save("manual_level2", "hemascribe_fine_sortgated",
                              "confusion_manual_level2_vs_hemascribe_fine")
    crosstab_save("population", "hemascribe_broad",
                  "confusion_sort_vs_hemascribe_broad_raw", normalize=None)
    crosstab_save("population", "hemascribe_fine",
                  "confusion_sort_vs_hemascribe_fine_raw", normalize=None)

    # ── Heatmap (manual_level2 × fine, top 20 fine labels by abundance) ─
    top_cols = (
        adata.obs["hemascribe_fine_sortgated"].value_counts().head(20).index.tolist()
    )
    sub = conf_fine.reindex(columns=[c for c in top_cols if c in conf_fine.columns])
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(sub.columns)),
                                    0.45 * len(sub.index) + 2.5))
    im = ax.imshow(sub.values, cmap="viridis", aspect="auto",
                   vmin=0, vmax=min(1.0, sub.values.max()))
    ax.set_xticks(range(len(sub.columns)))
    ax.set_xticklabels(sub.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(sub.index)))
    ax.set_yticklabels(sub.index)
    ax.set_xlabel("HemaScribe fine label (sort-gated)")
    ax.set_ylabel("manual_level2")
    ax.set_title("Row-normalised confusion (top 20 HemaScribe fine labels)")
    plt.colorbar(im, ax=ax, label="fraction of row")
    plt.tight_layout()
    fig.savefig(os.path.join(args.audit_dir, "confusion_heatmap_level2.png"),
                dpi=150)
    plt.close(fig)

    # ── README ──────────────────────────────────────────────────────────
    with open(os.path.join(args.audit_dir, "README.md"), "w") as fh:
        fh.write("# HemaScribe annotation audit\n\n")
        fh.write(f"Cells: {n}\n\n")
        fh.write(f"## Sort-gate conflict rates\n")
        fh.write(f"- broad: {n_flag_broad}/{n} ({n_flag_broad/n*100:.1f}%)\n")
        fh.write(f"- fine : {n_flag_fine}/{n} ({n_flag_fine/n*100:.1f}%)\n")
        fh.write(f"- any  : {n_flag_any}/{n} ({n_flag_any/n*100:.1f}%)\n\n")
        fh.write("## Files\n")
        fh.write("- confusion_manual_level1_vs_hemascribe_broad.csv (row-normalised)\n")
        fh.write("- confusion_manual_level2_vs_hemascribe_fine.csv (row-normalised)\n")
        fh.write("- confusion_sort_vs_hemascribe_broad_raw.csv (counts)\n")
        fh.write("- confusion_sort_vs_hemascribe_fine_raw.csv (counts)\n")
        fh.write("- sortgate_flag_rates.csv\n")
        fh.write("- confusion_heatmap_level2.png\n\n")
        fh.write("## Output AnnData columns\n")
        fh.write("- hemascribe_broad, hemascribe_fine, hemascribe_score (raw)\n")
        fh.write("- hemascribe_broad_sortgated, hemascribe_fine_sortgated (NaN where prior violated)\n")
        fh.write("- hemascribe_sortgate_flag (True where any granularity violates sort prior)\n")

    adata.write_h5ad(args.out)
    print(f"\n[21c] Wrote {args.out}")
    print(f"[21c] Audit in {args.audit_dir}/")


if __name__ == "__main__":
    main()
