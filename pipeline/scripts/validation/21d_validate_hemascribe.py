#!/usr/bin/env python3
"""
21d_validate_hemascribe.py

Three orthogonal sanity checks for the HemaScribe labels:

  (1) Marker-gene expression per HemaScribe label — does each label show the
      canonical mouse-BM transcriptional signature?
  (2) hematopoietic.score distribution per broad label — low-score calls are
      lower-confidence; flag biologically incongruent score/label pairs.
  (3) Triple agreement manual_level1 / popv_prediction / hemascribe_broad
      mapped to a coarse compartment vocabulary.

Specifically targets the Basophil_prog → fine_GMP cross-mapping concern:
the per-label marker dotplot for fine GMP shows whether those cells carry
GMP markers (Elane/Mpo/Prtn3) or basophil markers (Mcpt8/Prss34/Cpa3).

Outputs (in --out):
  marker_dotplot_broad.png
  marker_dotplot_fine.png
  marker_dotplot_baso_vs_gmp_focus.png   per-cell marker check on fine_GMP
  score_distribution_per_broad.png
  score_summary_per_broad.csv
  triple_agreement_coarse.csv
  triple_agreement_summary.md
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


# ── Canonical mouse BM markers (one panel; scanpy handles missing genes) ────
MARKERS = {
    "HSC":          ["Procr", "Hlf", "Mecom", "Hoxb5", "Mpl", "Hoxa9"],
    "MPP/early":    ["Cd34", "Kit", "Flt3", "Ly6a", "Slamf1"],
    "GMP":          ["Elane", "Mpo", "Prtn3", "Cebpa", "Cebpe", "Csf1r"],
    "Basophil":     ["Mcpt8", "Prss34", "Cpa3", "Cd200r3", "Lmo4", "Hgf"],
    "Erythroid":    ["Klf1", "Gata1", "Hbb-bt", "Hba-a1", "Gypa"],
    "Megakaryo":    ["Pf4", "Itga2b", "Vwf", "Gp1bb"],
    "Neutrophil":   ["S100a8", "S100a9", "Ltf", "Mmp9", "Camp"],
    "Monocyte":     ["Ly6c2", "Cd14", "Csf1r", "Ccr2"],
    "DC":           ["Siglech", "Bst2", "Itgax", "Irf8"],
    "B":            ["Cd19", "Ms4a1", "Ighm", "Pax5"],
    "T/NK":         ["Cd3d", "Cd3e", "Klrb1c", "Ncr1"],
    "Lymph_prog":   ["Il7r", "Dntt", "Vpreb1"],
}
MARKER_PANEL = [g for panel in MARKERS.values() for g in panel]


# ── Coarse compartment map (collapses each method's vocabulary) ─────────────
COARSE = {
    # stem / multipotent
    **dict.fromkeys(["LTHSC", "HSC", "STHSC", "MPP", "MPP1", "MPP2", "MPP3",
                     "MPP4", "LMPP", "HSPC", "FcG_neg_MPP3", "FcG_pos_MPP3",
                     "progenitor"], "Stem_MPP"),
    # myeloid progenitor
    **dict.fromkeys(["CMP", "GMP", "cMoP", "mGMP", "GP",
                     "myeloid_progenitor", "ambiguous_myeloid",
                     "cycling_myeloid"], "Myeloid_prog"),
    # erythroid / mega
    **dict.fromkeys(["MEP", "EryP", "MkP", "Megakaryocyte",
                     "erythroid"], "Ery_Mk"),
    # mature myeloid
    **dict.fromkeys(["Neutrophil", "Immature_neutrophil", "Monocyte",
                     "Mast_cell", "Basophil", "cDC", "pDC",
                     "granulocyte_neutrophil", "monocyte", "basophil",
                     "DC"], "Mature_myeloid"),
    # lymphoid
    **dict.fromkeys(["B_cell", "Immature_B_cell", "T_cell", "NK_cell",
                     "Plasma_cell", "CLP", "plasma_cell"], "Lymphoid"),
    # not-hematopoietic
    **dict.fromkeys(["NotHem", "RBC"], "NotHem"),
}


def to_coarse(s):
    return s.astype(str).map(lambda x: COARSE.get(x, "Other"))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default="results/21_hemascribe/validation")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"[21d] Loading {args.input}")
    adata = sc.read_h5ad(args.input)

    # Use lognorm layer for expression-based plots so genes are comparable.
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"]

    panel = [g for g in MARKER_PANEL if g in adata.var_names]
    missing = sorted(set(MARKER_PANEL) - set(panel))
    if missing:
        print(f"      missing from var: {missing}")
    print(f"      using {len(panel)} marker genes")

    # ── (1) Marker dotplot per HemaScribe broad ─────────────────────────
    print("[21d] Marker dotplot — hemascribe_broad")
    sc.pl.dotplot(
        adata, var_names={k: [g for g in v if g in adata.var_names]
                          for k, v in MARKERS.items()},
        groupby="hemascribe_broad", standard_scale="var",
        show=False, swap_axes=False,
    )
    plt.savefig(os.path.join(args.out, "marker_dotplot_broad.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    print("[21d] Marker dotplot — hemascribe_fine")
    sc.pl.dotplot(
        adata, var_names={k: [g for g in v if g in adata.var_names]
                          for k, v in MARKERS.items()},
        groupby="hemascribe_fine", standard_scale="var",
        show=False,
    )
    plt.savefig(os.path.join(args.out, "marker_dotplot_fine.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    # ── Focused diagnostic: fine_GMP cells split by manual_level2 ───────
    print("[21d] Focused GMP/Basophil marker check on fine_GMP cells")
    gmp_cells = adata[adata.obs["hemascribe_fine"] == "GMP"].copy()
    gmp_cells.obs["origin_manual_level2"] = (
        gmp_cells.obs["manual_level2"].astype(str)
    )
    if gmp_cells.n_obs > 0:
        focus_markers = (MARKERS["GMP"] + MARKERS["Basophil"]
                         + MARKERS["MPP/early"][:3])
        focus_markers = [g for g in focus_markers if g in adata.var_names]
        sc.pl.dotplot(
            gmp_cells, var_names=focus_markers,
            groupby="origin_manual_level2", standard_scale="var", show=False,
        )
        plt.savefig(
            os.path.join(args.out, "marker_dotplot_baso_vs_gmp_focus.png"),
            dpi=150, bbox_inches="tight",
        )
        plt.close("all")

    # ── (2) hematopoietic.score distribution per broad label ────────────
    print("[21d] Score distribution per hemascribe_broad")
    sub = adata.obs[["hemascribe_broad", "hemascribe_score"]].dropna()
    order = sub.groupby("hemascribe_broad", observed=True)["hemascribe_score"]\
              .median().sort_values(ascending=False).index.tolist()
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(order)), 5))
    data = [sub.loc[sub["hemascribe_broad"] == lab, "hemascribe_score"].values
            for lab in order]
    ax.boxplot(data, tick_labels=order, showfliers=False)
    ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel("hematopoietic.score")
    ax.set_title("HemaScribe hematopoietic.score per broad label")
    ax.axhline(0.10, color="red", linestyle="--", lw=0.8,
               label="low-confidence threshold")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "score_distribution_per_broad.png"),
                dpi=150)
    plt.close(fig)

    summary = sub.groupby("hemascribe_broad", observed=True)["hemascribe_score"]\
                 .describe().round(3)
    summary["pct_low_conf"] = (
        sub.assign(low=lambda d: d["hemascribe_score"] < 0.10)
           .groupby("hemascribe_broad", observed=True)["low"].mean()
           .round(3) * 100
    )
    summary.to_csv(os.path.join(args.out, "score_summary_per_broad.csv"))
    print(summary.to_string())

    # ── (3) Triple agreement on coarse compartments ─────────────────────
    print("[21d] Triple-agreement: manual / popV / HemaScribe (coarse)")
    obs = adata.obs.copy()
    obs["coarse_manual"]     = to_coarse(obs["manual_level1"])
    obs["coarse_popv"]       = to_coarse(obs["popv_prediction"])
    obs["coarse_hemascribe"] = to_coarse(obs["hemascribe_broad"])

    agree_3 = (
        (obs["coarse_manual"] == obs["coarse_popv"]) &
        (obs["coarse_popv"] == obs["coarse_hemascribe"])
    )
    agree_man_hema = obs["coarse_manual"] == obs["coarse_hemascribe"]
    agree_man_pop  = obs["coarse_manual"] == obs["coarse_popv"]
    agree_pop_hema = obs["coarse_popv"]   == obs["coarse_hemascribe"]

    n = len(obs)
    rows = [
        ("manual ∩ popV ∩ HemaScribe", agree_3.sum()),
        ("manual = HemaScribe",        agree_man_hema.sum()),
        ("manual = popV",              agree_man_pop.sum()),
        ("popV   = HemaScribe",        agree_pop_hema.sum()),
    ]
    tbl = pd.DataFrame(rows, columns=["pair", "n_agree"])
    tbl["pct"] = (tbl["n_agree"] / n * 100).round(1)
    tbl.to_csv(os.path.join(args.out, "triple_agreement_coarse.csv"),
               index=False)
    print(tbl.to_string(index=False))

    # Per-coarse-cluster breakdown
    per_compartment = (
        obs.groupby("coarse_manual", observed=True).apply(
            lambda d: pd.Series({
                "n_cells":           len(d),
                "agree_3way_pct":    100 * ((d["coarse_manual"] == d["coarse_popv"]) &
                                             (d["coarse_popv"] == d["coarse_hemascribe"])).mean(),
                "agree_manual_pop":  100 * (d["coarse_manual"] == d["coarse_popv"]).mean(),
                "agree_manual_hema": 100 * (d["coarse_manual"] == d["coarse_hemascribe"]).mean(),
                "agree_pop_hema":    100 * (d["coarse_popv"]   == d["coarse_hemascribe"]).mean(),
            })
        ).round(1).sort_values("n_cells", ascending=False)
    )
    per_compartment.to_csv(os.path.join(args.out,
                                        "triple_agreement_per_compartment.csv"))
    print("\nPer-compartment agreement:")
    print(per_compartment.to_string())

    with open(os.path.join(args.out, "triple_agreement_summary.md"), "w") as fh:
        fh.write(f"# Triple-agreement (coarse compartments)\n\n")
        fh.write(f"Total cells: {n}\n\n")
        fh.write("```\n" + tbl.to_string(index=False) + "\n```\n\n")
        fh.write("## Per-compartment\n\n")
        fh.write("```\n" + per_compartment.to_string() + "\n```\n")

    print(f"\n[21d] Validation outputs in {args.out}/")


if __name__ == "__main__":
    main()
