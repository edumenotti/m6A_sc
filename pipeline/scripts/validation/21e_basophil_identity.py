#!/usr/bin/env python3
"""
21e_basophil_identity.py

Investigate basophil-lineage labels: are HemaScribe's "Basophil" calls in the LK
gate (252 cells) and our manual "Basophil_prog" cluster (711 cells) truly basophil
PROGENITORS (BMCP, Cpa3+/Lmo4+ but Mcpt8-/Prss34-) or mature basophils that
leaked through the sort?

Per-cell marker logic:
- Mature basophil:   Mcpt8+/Prss34+/Cpa3+/Cd200r3+/Lmo4+   AND  Mpo-/Elane-/Cd34-
- BMCP / Baso_prog:  Cpa3+/Lmo4+/Hgf+                      AND  Mpo+/Cd34+  (transition state)
- GMP-like:          Mpo+/Elane+/Prtn3+/Cebpa+              AND  Cpa3-/Lmo4-

Outputs (in --out):
  baso_marker_scores_per_subgroup.csv   — mean expression + %+ per subgroup
  baso_dotplot_subgroups.png            — dotplot of canonical markers
  baso_classification_per_cell.csv      — per-cell verdict (mature / BMCP / GMP-like / ambiguous)
  baso_verdict_summary.md               — narrative
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


MARKERS = {
    "Mature_baso":    ["Mcpt8", "Prss34", "Cd200r3"],
    "BMCP_shared":    ["Cpa3", "Lmo4", "Hgf"],
    "GMP_program":    ["Mpo", "Elane", "Prtn3", "Cebpa"],
    "Stemness":       ["Kit", "Cd34", "Ly6a"],
}
ALL = [g for v in MARKERS.values() for g in v]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default="results/21_hemascribe/validation/basophil")
    return p.parse_args()


def classify(row):
    """Decide whether a cell is mature basophil, BMCP, GMP-like, or ambiguous."""
    mature = row["Mature_baso_pos"] >= 2
    bmcp   = row["BMCP_shared_pos"] >= 2
    gmp    = row["GMP_program_pos"] >= 2
    stem   = row["Stemness_pos"] >= 1

    if mature and not gmp:
        return "mature_basophil"
    if bmcp and gmp:
        return "BMCP_transition"
    if bmcp and not gmp:
        return "early_basophil"
    if gmp and not bmcp:
        return "GMP_like"
    return "ambiguous"


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"[21e] Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"]

    available = [g for g in ALL if g in adata.var_names]
    print(f"      using {len(available)} markers (missing: {sorted(set(ALL) - set(available))})")

    # Define the cell subgroups we want to compare
    sub_labels = pd.Series("other", index=adata.obs_names, dtype=object)
    # Manual Basophil_prog cluster (entire)
    sub_labels[adata.obs["manual_level2"] == "Basophil_prog"] = "manual_Basophil_prog"
    # Manual mature Basophil cluster
    sub_labels[adata.obs["manual_level2"] == "Basophil"]      = "manual_Basophil_mature"
    # HemaScribe Basophil call within LK gate (the 252)
    mask = (adata.obs["hemascribe_broad"] == "Basophil") & (adata.obs["population"] == "LK")
    sub_labels[mask] = "hemascribe_Baso_in_LK"
    # HemaScribe Basophil call within I gate (the 43 — mature)
    mask = (adata.obs["hemascribe_broad"] == "Basophil") & (adata.obs["population"] == "I")
    sub_labels[mask] = "hemascribe_Baso_in_I"
    # Reference: GMP_neutrophil_primed cells for comparison
    sub_labels[adata.obs["manual_level2"] == "GMP_neutrophil_primed"] = "ref_GMP_primed"
    adata.obs["baso_subgroup"] = pd.Categorical(sub_labels)

    targets = ["manual_Basophil_prog", "manual_Basophil_mature",
               "hemascribe_Baso_in_LK", "hemascribe_Baso_in_I",
               "ref_GMP_primed"]
    print("\n[21e] subgroup sizes:")
    print(adata.obs["baso_subgroup"].value_counts().reindex(targets).to_string())

    # ── Per-cell positivity per marker (>0 in lognorm = expressed) ──
    X = adata[:, available].X
    expr = pd.DataFrame((X > 0).toarray(), index=adata.obs_names, columns=available)
    counts = {}
    for grp, genes in MARKERS.items():
        gs = [g for g in genes if g in available]
        counts[f"{grp}_pos"] = expr[gs].sum(axis=1)
    counts_df = pd.DataFrame(counts)
    adata.obs = pd.concat([adata.obs, counts_df], axis=1)

    # ── Per-cell classification ──
    verdict = counts_df.apply(classify, axis=1)
    adata.obs["baso_verdict"] = pd.Categorical(verdict)

    # ── Per-subgroup mean/percent table ──
    rows = []
    for grp in targets:
        cells = adata.obs.index[adata.obs["baso_subgroup"] == grp]
        if len(cells) == 0: continue
        row = {"subgroup": grp, "n_cells": len(cells)}
        for gene in available:
            sub_X = adata[cells, gene].X.toarray().ravel()
            row[f"{gene}_mean"] = float(np.mean(sub_X))
            row[f"{gene}_pct_pos"] = float((sub_X > 0).mean() * 100)
        for col in counts_df.columns:
            row[col + "_mean"] = float(counts_df.loc[cells, col].mean())
        # verdict distribution
        v = verdict.loc[cells].value_counts(normalize=True).round(3).to_dict()
        for k in ["mature_basophil", "BMCP_transition", "early_basophil",
                  "GMP_like", "ambiguous"]:
            row[f"verdict_{k}"] = float(v.get(k, 0.0))
        rows.append(row)
    out_table = pd.DataFrame(rows).set_index("subgroup")
    out_table.to_csv(os.path.join(args.out, "baso_marker_scores_per_subgroup.csv"))

    # Compact verdict table
    verdict_cols = [c for c in out_table.columns if c.startswith("verdict_")]
    pct_cols     = [c for c in out_table.columns if c.endswith("_pct_pos")]
    print("\n[21e] Verdict distribution per subgroup:")
    print(out_table[verdict_cols].round(2).to_string())
    print("\n[21e] %-positive per marker (per subgroup):")
    print(out_table[pct_cols].round(0).to_string())

    # ── Dotplot ──
    sc.pl.dotplot(adata, var_names={k: [g for g in v if g in available]
                                     for k, v in MARKERS.items()},
                  groupby="baso_subgroup", standard_scale="var", show=False)
    plt.savefig(os.path.join(args.out, "baso_dotplot_subgroups.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    # ── Per-cell verdict csv ──
    pcd = pd.DataFrame({
        "subgroup":     adata.obs["baso_subgroup"],
        "verdict":      adata.obs["baso_verdict"],
        "population":   adata.obs["population"],
        "manual_level2":adata.obs["manual_level2"],
        "hemascribe_broad": adata.obs["hemascribe_broad"],
    })
    pcd[pcd["subgroup"] != "other"].to_csv(
        os.path.join(args.out, "baso_classification_per_cell.csv"))

    # ── Narrative ──
    with open(os.path.join(args.out, "baso_verdict_summary.md"), "w") as fh:
        fh.write("# Basophil lineage validation\n\n")
        fh.write("## Verdict distribution per subgroup\n\n```\n")
        fh.write(out_table[verdict_cols].round(2).to_string())
        fh.write("\n```\n\n## %-positive per marker per subgroup\n\n```\n")
        fh.write(out_table[pct_cols].round(0).to_string())
        fh.write("\n```\n")

    print(f"\n[21e] Outputs in {args.out}/")


if __name__ == "__main__":
    main()
