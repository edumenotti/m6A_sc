#!/usr/bin/env python3
"""
20_pathway_activity.py

Pathway activity per celltype × contrast, computed from the PyDESeq2 result
tables in `results/19_pseudobulk_deg/`. We use decoupler's univariate linear
model (ULM) on the `stat` column (the Wald statistic) — recommended pattern
for transferring bulk-style DEG output into a multi-collection enrichment.

Collections (mouse):
- MSigDB Hallmark
- MSigDB Reactome pathways
- MSigDB KEGG pathways

Pathway scoring with n=2 is more robust than individual gene calls because
the test aggregates ~50–500 genes per set.

Inputs
------
--deg-dir     directory with deg_<celltype>_<contrast>.csv (from script 19)
--out         output directory
--organism    'mouse' (default) or 'human'
--padj        threshold for "significant" pathway hits (default 0.05)

Outputs
-------
pathway_activity_<contrast>.csv         long-format scores/padj (all celltypes)
pathway_heatmap_<contrast>.png          celltypes × pathway score heatmap
pathway_top_hits.csv                    union of significant hits across runs
"""
import argparse
import glob
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--deg-dir", required=True)
    p.add_argument("--out", default="results/20_pathway_activity")
    p.add_argument("--organism", default="mouse", choices=["mouse", "human"])
    p.add_argument("--padj", type=float, default=0.05)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    import decoupler as dc

    # ── Build prior knowledge network from MSigDB ───────────────────────
    # NOTE: decoupler 2.1.6 exposes MSigDB via `dc.op.resource('MSigDB', ...)`.
    # The returned DataFrame has columns ['genesymbol','collection','geneset'].
    # We filter the single resource by the `collection` field to assemble
    # Hallmark + Reactome + KEGG, then rename to source/target/weight.
    print("Fetching MSigDB collections ...")
    msig = dc.op.resource("MSigDB", organism=args.organism)
    keep_collections = {
        "hallmark",
        "reactome_pathways",
        "kegg_pathways",
        "kegg_medicus_pathways",
    }
    net = msig[msig["collection"].isin(keep_collections)].copy()
    net = net.rename(columns={"geneset": "source", "genesymbol": "target"})
    net["weight"] = 1.0
    net = net[["source", "target", "weight"]].drop_duplicates()
    print(
        f"  Combined network: {net['source'].nunique()} pathways, "
        f"{len(net)} edges"
    )

    # ── Load DEG result tables ──────────────────────────────────────────
    files = sorted(glob.glob(os.path.join(args.deg_dir, "deg_*_treatment.csv"))) + \
            sorted(glob.glob(os.path.join(args.deg_dir, "deg_*_genotype.csv")))
    files = [f for f in files if "_top" not in f]
    if not files:
        raise SystemExit(f"No DEG CSVs found in {args.deg_dir}")

    # Long-format: celltype × pathway scores per contrast
    all_results = []
    for f in files:
        base = os.path.basename(f).replace("deg_", "").replace(".csv", "")
        parts = base.rsplit("_", 1)
        celltype, contrast = parts[0], parts[1]
        df = pd.read_csv(f).dropna(subset=["stat"])
        if df.empty:
            continue
        mat = df.set_index("gene")[["stat"]].T   # 1 × n_genes
        try:
            # dc.mt.ulm returns (activities, RAW p-values) — not FDR-corrected.
            acts, pvals = dc.mt.ulm(
                data=mat, net=net, tmin=5, verbose=False
            )
        except Exception as e:
            print(f"  [{celltype}/{contrast}] ULM failed: {e}")
            continue
        long = pd.DataFrame({
            "pathway": acts.columns,
            "score":   acts.values.ravel(),
            "pvalue":  pvals.values.ravel(),
        })
        long["celltype"] = celltype
        long["contrast"] = contrast
        all_results.append(long)

    if not all_results:
        raise SystemExit("No pathway results produced — check DEG inputs")

    big = pd.concat(all_results, ignore_index=True)

    # Apply BH FDR per (celltype, contrast) — each (celltype, contrast) tests
    # ~1500 pathways. Correcting within that family is the standard decoupler
    # pattern and avoids the false-positive inflation flagged in review.
    from statsmodels.stats.multitest import multipletests
    big["padj"] = np.nan
    for (ct, cn), idx in big.groupby(["celltype", "contrast"]).groups.items():
        mask = big.loc[idx, "pvalue"].notna()
        if not mask.any():
            continue
        pv = big.loc[idx[mask], "pvalue"].values
        _, padj_vals, _, _ = multipletests(pv, method="fdr_bh")
        big.loc[idx[mask], "padj"] = padj_vals

    # ── Per-contrast wide tables + heatmaps ─────────────────────────────
    for contrast in big["contrast"].unique():
        sub = big[big["contrast"] == contrast].copy()
        sub.to_csv(
            os.path.join(args.out, f"pathway_activity_{contrast}.csv"),
            index=False
        )

        # Pick top pathways by max |score| across celltypes
        ranking = sub.groupby("pathway")["score"].apply(
            lambda x: x.abs().max()
        ).sort_values(ascending=False)
        top_paths = ranking.head(40).index.tolist()
        heat = sub[sub["pathway"].isin(top_paths)].pivot_table(
            index="pathway", columns="celltype", values="score", fill_value=0
        ).reindex(top_paths)

        fig, ax = plt.subplots(figsize=(max(8, 0.6*heat.shape[1]),
                                        max(8, 0.3*heat.shape[0])))
        vmax = float(np.nanmax(np.abs(heat.values))) or 1.0
        im = ax.imshow(heat.values, cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(heat.shape[1]))
        ax.set_xticklabels(heat.columns, rotation=45, ha="right")
        ax.set_yticks(range(heat.shape[0]))
        ax.set_yticklabels(heat.index, fontsize=7)
        ax.set_title(f"Top 40 pathway activities — contrast: {contrast}")
        plt.colorbar(im, ax=ax, label="ULM score")
        plt.tight_layout()
        fig.savefig(
            os.path.join(args.out, f"pathway_heatmap_{contrast}.png"),
            dpi=150
        )
        plt.close(fig)

    # ── Union of significant pathway hits ───────────────────────────────
    sig = big[big["padj"] < args.padj].copy()
    sig = sig.sort_values(["contrast", "celltype", "padj"])
    sig.to_csv(os.path.join(args.out, "pathway_top_hits.csv"), index=False)

    print(f"\nScript 20 complete. Outputs in {args.out}/")


if __name__ == "__main__":
    main()
