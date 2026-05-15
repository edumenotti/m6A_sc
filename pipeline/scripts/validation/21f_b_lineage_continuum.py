#!/usr/bin/env python3
"""
21f_b_lineage_continuum.py

Validate HemaScribe's split of manual Pro_Pre_B (343 cells) into B_cell (50%) +
Immature_B_cell (38%) + CLP (4%). Is this a real maturation continuum or noise?

Canonical mouse B-lineage maturation:
  CLP       → ProB         → PreB          → ImmatureB   → MatureB
  Il7r+     Cd34+/Rag1+    Vpreb1+/Igll1+   IgM+         IgM+/Cd19+
  Flt3+     Vpreb1+/Igll1+ Cd19+/Rag1+      Cd19+        Ms4a1+
  Cd34+     Cd19-

Markers:
- CLP:           Il7r, Flt3, Cd34, Dntt
- ProB:          Rag1, Vpreb1, Igll1, Cd34 (Cd19-)
- PreB:          Vpreb1, Igll1, Cd19, Rag1
- ImmatureB:     Ighm, Cd19, Ms4a1, Cd24a
- MatureB:       Cd19, Ms4a1, Ighd, Cd23 (Fcer2a)

Outputs:
  b_lineage_dotplot.png
  b_lineage_marker_table.csv
  b_lineage_verdict.md
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
    "CLP":      ["Il7r", "Flt3", "Cd34", "Dntt"],
    "ProB":     ["Rag1", "Vpreb1", "Igll1"],
    "PreB":     ["Vpreb1", "Igll1", "Cd19", "Rag1"],
    "Immature": ["Ighm", "Cd19", "Ms4a1", "Cd24a"],
    "Mature":   ["Cd19", "Ms4a1", "Ighd", "Fcer2a", "Pax5"],
}
ALL = sorted({g for v in MARKERS.values() for g in v})


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default="results/21_hemascribe/validation/b_lineage")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"[21f] Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"]

    available = [g for g in ALL if g in adata.var_names]
    print(f"      markers found: {len(available)}/{len(ALL)} "
          f"(missing: {sorted(set(ALL) - set(available))})")

    # Subgroups: split manual Pro_Pre_B by HemaScribe broad call
    obs = adata.obs.copy()
    sub = pd.Series("other", index=obs.index, dtype=object)

    ppb_mask = obs["manual_level2"] == "Pro_Pre_B"
    sub[ppb_mask & (obs["hemascribe_broad"] == "B_cell")]         = "ProPreB→Bcell"
    sub[ppb_mask & (obs["hemascribe_broad"] == "Immature_B_cell")] = "ProPreB→ImmB"
    sub[ppb_mask & (obs["hemascribe_broad"] == "CLP")]              = "ProPreB→CLP"
    sub[ppb_mask & ~obs["hemascribe_broad"].isin(
        ["B_cell", "Immature_B_cell", "CLP"])]                       = "ProPreB→other"
    # Reference groups
    sub[obs["manual_level2"] == "Mature_B"] = "ref_Mature_B"
    adata.obs["b_subgroup"] = pd.Categorical(sub)

    targets = ["ProPreB→CLP", "ProPreB→ImmB", "ProPreB→Bcell",
               "ProPreB→other", "ref_Mature_B"]
    print("\n[21f] subgroup sizes:")
    print(adata.obs["b_subgroup"].value_counts().reindex(targets, fill_value=0).to_string())

    # Marker mean + %-positive per subgroup
    rows = []
    for grp in targets:
        cells = adata.obs.index[adata.obs["b_subgroup"] == grp]
        if len(cells) == 0: continue
        row = {"subgroup": grp, "n_cells": len(cells)}
        for g in available:
            x = adata[cells, g].X.toarray().ravel()
            row[f"{g}_mean"]    = float(x.mean())
            row[f"{g}_pct_pos"] = float((x > 0).mean() * 100)
        rows.append(row)
    table = pd.DataFrame(rows).set_index("subgroup")
    table.to_csv(os.path.join(args.out, "b_lineage_marker_table.csv"))

    # Console summary: %-positive only
    pct_cols = [c for c in table.columns if c.endswith("_pct_pos")]
    print("\n[21f] %-positive per marker:")
    print(table[pct_cols].round(0).to_string())

    sc.pl.dotplot(adata, var_names={k: [g for g in v if g in available]
                                     for k, v in MARKERS.items()},
                  groupby="b_subgroup", standard_scale="var", show=False)
    plt.savefig(os.path.join(args.out, "b_lineage_dotplot.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    with open(os.path.join(args.out, "b_lineage_verdict.md"), "w") as fh:
        fh.write("# B-lineage continuum validation\n\n")
        fh.write("Manual Pro_Pre_B was split by HemaScribe into B_cell + Immature_B + CLP.\n")
        fh.write("If real maturation: CLP→ProPreB→ImmB→Bcell should show ordered "
                 "marker progression (Dntt/Flt3↓, Rag1/Vpreb1 mid, Cd19/Ms4a1↑, Ighd late).\n\n")
        fh.write("## %-positive per marker per subgroup\n\n```\n")
        fh.write(table[pct_cols].round(0).to_string())
        fh.write("\n```\n")
    print(f"\n[21f] Outputs in {args.out}/")


if __name__ == "__main__":
    main()
