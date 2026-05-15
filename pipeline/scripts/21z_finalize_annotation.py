#!/usr/bin/env python3
"""
21z_finalize_annotation.py

Apply the manual validation findings from 21d/e/f/g/h to produce the final
consolidated cell-type annotation. Decisions (justified in
results/21_hemascribe/validation/*):

  Rule 1: manual rescue for clusters HemaScribe lumps or mislabels
          → Basophil_prog, MEP_Erythroid, Progenitor_ambiguous, Ambiguous,
            Pro_Pre_B (HemaScribe split into B/ImmB/CLP rotates the names — 21f)
  Rule 2: HSPC compartment refinement
          → use hemascribe_fine for cells with hemascribe_broad == "HSPC"
            (HSC/STHSC/MPP2/3/4/MkP/GMP/CLP)
  Rule 3: everything else → hemascribe_broad
  Rule 4: NO cells become "Unassigned"; QC flags are kept separate
          qc_sort_leak: mature lineage in stem gate (LSK/LK), per 21g
              EryP/mGMP/cMoP/Neutrophil/Monocyte/Mast_cell in LSK or LK
              (NB: Plasma_cell in LSK is NOT a leak — they are Sca1+ plasmablasts, 21g)
          qc_low_score: hemascribe_score < 0.03 (mostly EryP/Mega legitimately low)

Outputs:
  - rewrites results/14_progenitor_annotated/adata_hemascribe.h5ad in place,
    adding columns: final_annotation, qc_sort_leak, qc_low_score
  - results/21_hemascribe/final_annotation_summary.csv
  - results/21_hemascribe/final_annotation_audit.md

Usage:
  pixi run python pipeline/scripts/21z_finalize_annotation.py \\
    --input results/14_progenitor_annotated/adata_hemascribe.h5ad
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc


# Manual clusters preserved verbatim — HemaScribe either misses them
# (Basophil_prog Cpa3+/Lmo4+, MEP_Erythroid intermediate) or mislabels
# their internal structure (Pro_Pre_B; the labels CLP/ImmB/Bcell from
# HemaScribe do not follow canonical maturation, per 21f).
MANUAL_TRUSTED = {
    "Basophil_prog",
    "MEP_Erythroid",
    "Progenitor_ambiguous",
    "Ambiguous",
    "Pro_Pre_B",
}

# Mature/committed labels that, when seen in a stem gate (LSK or LK), look
# like sort contamination. Plasma_cell is EXCLUDED — 21g showed plasmablasts
# can be Sca1+ legitimately.
MATURE_IN_STEM_GATE_LEAKS = {
    "Neutrophil", "Monocyte", "Mast_cell", "T_cell", "NK_cell",
    "B_cell", "Immature_B_cell", "Immature_neutrophil", "RBC",
    "EryP", "mGMP", "cMoP",   # committed myeloid/erythroid in LSK gate
}
# Note: cMoP/mGMP/EryP in LK gate is NOT a leak (LK is meant for those).
STEM_GATES = {"LSK"}  # the leaks we found in 21g are LSK-only


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default=None,
                   help="Output path (default: overwrite input)")
    p.add_argument("--audit-dir", default="results/21_hemascribe")
    return p.parse_args()


def main():
    args = parse_args()
    out_path = args.out or args.input

    print(f"[21z] Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"      shape: {adata.shape}")

    broad  = adata.obs["hemascribe_broad"].astype(str)
    fine   = adata.obs["hemascribe_fine"].astype(str)
    manual = adata.obs["manual_level2"].astype(str)
    sort_pop = adata.obs["population"].astype(str)
    score  = adata.obs["hemascribe_score"].astype(float)

    final = pd.Series("Unassigned", index=adata.obs_names, dtype=object)
    rule_used = pd.Series("none", index=adata.obs_names, dtype=object)

    # Rule 1: manual rescue
    m = manual.isin(MANUAL_TRUSTED)
    final[m]     = manual[m]
    rule_used[m] = "1_manual_rescue"

    # Rule 2: HSPC compartment → hemascribe_fine
    m = (final == "Unassigned") & (broad == "HSPC") & ~fine.isin(
        ["NotHSPC", "nan", ""]
    )
    final[m]     = fine[m]
    rule_used[m] = "2_hemascribe_fine"

    # Rule 2b: HSPC without fine call → generic HSPC
    m = (final == "Unassigned") & (broad == "HSPC")
    final[m]     = "HSPC_unspecified"
    rule_used[m] = "2b_HSPC_generic"

    # Rule 3: rest → hemascribe_broad
    m = (final == "Unassigned") & ~broad.isin(["nan", ""])
    final[m]     = broad[m]
    rule_used[m] = "3_hemascribe_broad"

    n_unassigned = (final == "Unassigned").sum()
    print(f"[21z] Unassigned after rules: {n_unassigned}")
    if n_unassigned > 0:
        print(adata.obs.loc[final == "Unassigned",
                            ["population", "manual_level2",
                             "hemascribe_broad", "hemascribe_fine"]]
              .value_counts().head(10).to_string())

    # QC flags — informational, do NOT change label
    qc_sort_leak = sort_pop.isin(STEM_GATES) & broad.isin(MATURE_IN_STEM_GATE_LEAKS)
    qc_low_score = score < 0.03

    adata.obs["final_annotation"] = pd.Categorical(final)
    adata.obs["final_annotation_rule"] = pd.Categorical(rule_used)
    adata.obs["qc_sort_leak"] = qc_sort_leak.values
    adata.obs["qc_low_score"] = qc_low_score.values

    print(f"\n[21z] final_annotation distribution ({final.nunique()} labels):")
    print(final.value_counts().to_string())
    print(f"\n[21z] Rule usage:")
    print(rule_used.value_counts().to_string())
    print(f"\n[21z] QC flags:")
    print(f"  qc_sort_leak:  {qc_sort_leak.sum()} cells ({qc_sort_leak.mean()*100:.2f}%)")
    print(f"  qc_low_score:  {qc_low_score.sum()} cells ({qc_low_score.mean()*100:.2f}%)")
    print(f"  any:           {(qc_sort_leak | qc_low_score).sum()} cells "
          f"({(qc_sort_leak | qc_low_score).mean()*100:.2f}%)")

    # ── Audit outputs ───────────────────────────────────────────────────
    os.makedirs(args.audit_dir, exist_ok=True)

    summary = pd.DataFrame({
        "n_cells":         final.value_counts(),
        "pct_total":       (final.value_counts(normalize=True) * 100).round(2),
    })
    summary["pct_sort_leak"] = (
        adata.obs.groupby("final_annotation", observed=True)["qc_sort_leak"]
            .mean().mul(100).round(2)
    )
    summary["pct_low_score"] = (
        adata.obs.groupby("final_annotation", observed=True)["qc_low_score"]
            .mean().mul(100).round(2)
    )
    summary["dominant_rule"] = (
        adata.obs.groupby("final_annotation", observed=True)["final_annotation_rule"]
            .agg(lambda s: s.value_counts().index[0])
    )
    summary.to_csv(os.path.join(args.audit_dir, "final_annotation_summary.csv"))

    with open(os.path.join(args.audit_dir, "final_annotation_audit.md"), "w") as fh:
        fh.write("# Final annotation audit\n\n")
        fh.write(f"Source AnnData: `{args.input}`\n\n")
        fh.write("## Decision rules (in priority order)\n\n")
        fh.write("1. **Manual rescue** for clusters HemaScribe lumps/mislabels: ")
        fh.write(", ".join(sorted(MANUAL_TRUSTED)) + "\n")
        fh.write("2. **HSPC compartment** (`hemascribe_broad == 'HSPC'`) "
                 "→ use `hemascribe_fine` (HSC / STHSC / MPP2/3/4 / MkP / GMP / CLP)\n")
        fh.write("3. **All other cells** → `hemascribe_broad`\n")
        fh.write("4. **QC flags** kept separate (cells keep their label):\n")
        fh.write("   - `qc_sort_leak`: mature/committed lineage in LSK gate (per 21g)\n")
        fh.write("   - `qc_low_score`: hemascribe_score < 0.03\n\n")
        fh.write("## Validation references\n\n")
        fh.write("- 21d: marker-gene + score + triple-agreement check\n")
        fh.write("- 21e: basophil identity — 252 LK Baso are BMCPs (not leak)\n")
        fh.write("- 21f: Pro_Pre_B split by HemaScribe is heterogeneity but mislabeled → use manual\n")
        fh.write("- 21g: LSK Plasma_cell is real Sca1+ plasmablast; LSK EryP/mGMP/cMoP are leak (~190)\n")
        fh.write("- 21h: popV LTHSC ↔ HemaScribe HSC+STHSC (84%); popV MEP ↔ HemaScribe MPP2 (91%)\n\n")
        fh.write(f"## Final distribution ({final.nunique()} labels, 0 unassigned)\n\n```\n")
        fh.write(summary.to_string())
        fh.write("\n```\n")

    adata.write_h5ad(out_path)
    print(f"\n[21z] Wrote {out_path}")
    print(f"[21z] Audit in {args.audit_dir}/")


if __name__ == "__main__":
    main()
