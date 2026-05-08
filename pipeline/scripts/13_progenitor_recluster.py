#!/usr/bin/env python
"""Progenitor subset recluster — diagnostic analysis for human review.

Extracts leiden_r1.0 clusters 3, 6, 10 (progenitors) plus 7, 8 (cycling
neighbors) from adata_annotated_final.h5ad. Re-embeds in progenitor-only PCA
space (IEG/mito/ribo genes excluded from HVGs). Computes Leiden clusterings at
multiple resolutions, runs DE tests, scores lineage markers. Saves UMAPs,
dotplots, DE tables, and an annotation_template.tsv for human review.

Usage:
    pixi run python pipeline/scripts/13_progenitor_recluster.py \\
        --input results/12_subset_recluster/adata_annotated_final.h5ad \\
        --out results/13_progenitor_recluster
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 1

# ── Cluster IDs ──────────────────────────────────────────────────────────────
PROG_CLUSTERS     = ["3", "6", "10"]        # primary progenitor clusters
NEIGHBOR_CLUSTERS = ["7", "8"]              # cycling myeloid — neighbor context
CLUSTER_KEY       = "leiden_r1.0"

# ── Resolutions to explore ───────────────────────────────────────────────────
RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]

# ── Genes to exclude from HVG selection and DE ranking ───────────────────────
IEG_GENES = {
    "Fos", "Fosb", "Fosl1", "Fosl2",
    "Jun", "Junb", "Jund",
    "Egr1", "Egr2", "Egr3",
    "Zfp36", "Zfp36l1", "Zfp36l2",
    "Btg2", "Nr4a1", "Nr4a2", "Nr4a3",
    "Dusp1", "Dusp5", "Dusp6",
    "Klf2", "Klf4", "Klf6",
    "Ier2", "Ier3", "Ier5",
    "Atf3", "Hspa1a", "Hspa1b",
}

def is_problem_gene(gene: str) -> bool:
    return (
        gene.startswith("mt-")
        or gene.startswith("Rps")
        or gene.startswith("Rpl")
        or gene in IEG_GENES
    )

# ── Progenitor lineage markers ────────────────────────────────────────────────
PROGENITOR_MARKERS: dict[str, list[str]] = {
    "HSC":           ["Hlf", "Procr", "Mecom", "Ly6a", "Slamf1", "Mpl"],
    "MPP":           ["Kit", "Cd34", "Adgrg1", "Adgrl4"],
    "LMPP":          ["Flt3", "Dntt", "Il7r"],
    "CMP":           ["Sox4", "Itga2b", "Gata2", "Angpt1"],
    "MEP":           ["Gata1", "Klf1", "Gypa", "Car1"],
    "GMP":           ["Mpo", "Elane", "Ctsg", "Prtn3", "Cebpa"],
    "Basophil_prog": ["Alox15", "Ms4a2", "Prss34", "Mcpt8", "Hdc"],
    "Cycling":       ["Mki67", "Top2a", "Pcna", "Birc5"],
    "Stress_IEG":    ["Fos", "Jun", "Egr1"],
}

QC_MARKERS = ["mt-Co1", "n_counts", "n_genes", "pct_counts_mt"]


# ── Helpers ───────────────────────────────────────────────────────────────────

def present(adata: sc.AnnData, genes: list[str]) -> list[str]:
    return [g for g in genes if g in adata.var_names]


def score_markers(adata: sc.AnnData) -> sc.AnnData:
    for name, genes in PROGENITOR_MARKERS.items():
        ok = present(adata, genes)
        if ok:
            sc.tl.score_genes(adata, ok, score_name=f"score_{name}", use_raw=False)
    return adata


def clean_and_reembed(
    adata_sub: sc.AnnData,
    n_hvg: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
) -> sc.AnnData:
    """HVG selection excluding problem genes, then PCA → neighbors → UMAP."""
    if "lognorm" in adata_sub.layers:
        adata_sub.X = adata_sub.layers["lognorm"].copy()
    # subset genes to non-problem genes for HVG selection
    clean_mask = ~pd.Series(adata_sub.var_names).apply(is_problem_gene).values
    adata_clean_var = adata_sub[:, clean_mask].copy()
    sc.pp.highly_variable_genes(adata_clean_var, n_top_genes=n_hvg, subset=False)
    hvg_names = adata_clean_var.var_names[adata_clean_var.var["highly_variable"]].tolist()
    # embed using those HVGs on the full adata (so X_pca stored in adata_sub)
    adata_sub.var["hvg_clean"] = adata_sub.var_names.isin(hvg_names)
    sc.pp.pca(adata_sub, n_comps=n_pcs, use_highly_variable=False,
              mask_var="hvg_clean")
    sc.pp.neighbors(adata_sub, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata_sub)
    return adata_sub


def leiden_multi(adata_sub: sc.AnnData, resolutions: list[float]) -> sc.AnnData:
    for res in resolutions:
        sc.tl.leiden(adata_sub, resolution=res, key_added=f"leiden_sub_r{res}")
    return adata_sub


def de_vs_background(
    adata_bg: sc.AnnData,
    leiden_key: str,
    out: Path,
    tag: str,
) -> pd.DataFrame:
    """DE each leiden cluster vs the rest of adata_bg, excluding problem genes."""
    clean_mask = ~pd.Series(adata_bg.var_names).apply(is_problem_gene).values
    a = adata_bg[:, clean_mask].copy()
    sc.tl.rank_genes_groups(
        a, groupby=leiden_key, method="wilcoxon",
        use_raw=False, key_added="rank_genes",
    )
    clusters = a.obs[leiden_key].cat.categories.tolist()
    rows = []
    for c in clusters:
        names  = a.uns["rank_genes"]["names"][c][:30]
        scores = a.uns["rank_genes"]["scores"][c][:30]
        pvals  = a.uns["rank_genes"]["pvals_adj"][c][:30]
        for rank, (g, s, p) in enumerate(zip(names, scores, pvals)):
            rows.append({"cluster": c, "rank": rank + 1, "gene": g,
                         "score": round(float(s), 3), "pval_adj": float(p)})
    df = pd.DataFrame(rows)
    df.to_csv(out / f"de_{tag}.csv", index=False)
    return df


def markers_in_top30(de_df: pd.DataFrame) -> pd.DataFrame:
    """For each cluster, count how many progenitor markers appear in top 30 DE."""
    results = []
    for cluster, grp in de_df.groupby("cluster"):
        top30 = set(grp.head(30)["gene"].tolist())
        row = {"cluster": cluster}
        for lineage, genes in PROGENITOR_MARKERS.items():
            hits = [g for g in genes if g in top30]
            row[lineage] = ", ".join(hits) if hits else ""
        results.append(row)
    return pd.DataFrame(results)


def save_umap(adata_sub: sc.AnnData, color: str, out: Path, prefix: str) -> None:
    key = color.replace("-", "_").replace("/", "_").replace(".", "")
    in_obs = color in adata_sub.obs.columns
    in_var = color in adata_sub.var_names
    if not in_obs and not in_var:
        return
    is_cat = in_obs and pd.api.types.is_categorical_dtype(adata_sub.obs[color])
    sc.pl.umap(
        adata_sub, color=color, show=False, frameon=False,
        legend_loc="right margin" if is_cat else "best",
        legend_fontsize=7,
    )
    plt.savefig(out / f"{prefix}_umap_{key}.png", dpi=150, bbox_inches="tight")
    plt.close()


def save_dotplot(
    adata_sub: sc.AnnData, groupby: str, out: Path, tag: str,
) -> None:
    genes_flat = [g for genes in PROGENITOR_MARKERS.values()
                  for g in present(adata_sub, genes)]
    genes_flat = list(dict.fromkeys(genes_flat))  # deduplicate, preserve order
    if not genes_flat:
        return
    try:
        sc.pl.dotplot(
            adata_sub, var_names=genes_flat, groupby=groupby,
            show=False, figsize=(max(8, len(genes_flat) * 0.4), 4),
        )
        plt.savefig(out / f"dotplot_{tag}.png", dpi=150, bbox_inches="tight")
    except Exception:
        pass
    plt.close()


def build_annotation_template(
    adata_sub: sc.AnnData,
    resolutions: list[float],
    marker_hits: dict[str, pd.DataFrame],
    out: Path,
) -> None:
    """Write template TSV for human to fill in."""
    rows = []
    for res in resolutions:
        key = f"leiden_sub_r{res}"
        df_hits = marker_hits.get(key, pd.DataFrame())
        for c in sorted(adata_sub.obs[key].cat.categories, key=int):
            mask = adata_sub.obs[key] == c
            n = int(mask.sum())
            orig_clusters = (
                adata_sub.obs.loc[mask, CLUSTER_KEY].value_counts().to_dict()
            )
            orig_str = " ".join(f"r1.0_{k}:{v}" for k, v in sorted(orig_clusters.items()))
            # top marker hits summary
            if not df_hits.empty:
                clust_row = df_hits.loc[df_hits["cluster"] == c]
                marker_summary = "; ".join(
                    f"{col}={val}" for col in PROGENITOR_MARKERS
                    for val in [clust_row.iloc[0][col] if not clust_row.empty else ""]
                    if val
                ) if not clust_row.empty else ""
            else:
                marker_summary = ""
            rows.append({
                "leiden_resolution": res,
                "subcluster_id": c,
                "n_cells": n,
                "source_clusters_r1.0": orig_str,
                "marker_hits": marker_summary,
                "new_level1": "",
                "new_level2": "",
                "confidence": "",
                "notes": "",
            })
    df = pd.DataFrame(rows)
    df.to_csv(out / "annotation_template.tsv", sep="\t", index=False)
    print(f"\n✓ Annotation template written to {out}/annotation_template.tsv")
    print("  Fill in: new_level1, new_level2, confidence, notes")
    print("  Then copy to pipeline/config/progenitor_annotation_map.tsv")
    print("  Also note which leiden_resolution you chose at top of that file (comment line)")


# ── Main ──────────────────────────────────────────────────────────────────────

def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # ── 1. Load ──────────────────────────────────────────────────────────────
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"obs columns: {list(adata.obs.columns)}")

    # ── 2. Subset ────────────────────────────────────────────────────────────
    prog_mask = adata.obs[CLUSTER_KEY].astype(str).isin(PROG_CLUSTERS)
    neighbor_mask = adata.obs[CLUSTER_KEY].astype(str).isin(NEIGHBOR_CLUSTERS)

    adata_prog = adata[prog_mask].copy()
    adata_prog_n = adata[prog_mask | neighbor_mask].copy()
    print(f"\nProgenitor subset (clusters {PROG_CLUSTERS}): {adata_prog.n_obs} cells")
    print(f"Progenitor + neighbor subset: {adata_prog_n.n_obs} cells")
    print("Cluster breakdown:")
    print(adata_prog.obs[CLUSTER_KEY].value_counts().sort_index().to_string())

    # Save original UMAP for reference plots
    adata_prog.obsm["X_umap_original"] = adata_prog.obsm["X_umap"].copy()
    adata_prog_n.obsm["X_umap_original"] = adata_prog_n.obsm["X_umap"].copy()

    # ── 3. Re-embed progenitor subset ────────────────────────────────────────
    print("\nRe-embedding progenitor subset (excluding mito/ribo/IEG from HVGs)...")
    adata_prog = clean_and_reembed(adata_prog)

    # ── 4. Leiden clusterings ────────────────────────────────────────────────
    print(f"\nComputing Leiden at resolutions: {RESOLUTIONS}")
    adata_prog = leiden_multi(adata_prog, RESOLUTIONS)
    for res in RESOLUTIONS:
        key = f"leiden_sub_r{res}"
        n_clusters = adata_prog.obs[key].nunique()
        print(f"  r={res}: {n_clusters} clusters")

    # Save barcodes + leiden assignments for downstream apply script
    leiden_cols = [CLUSTER_KEY] + [f"leiden_sub_r{r}" for r in RESOLUTIONS]
    adata_prog.obs[leiden_cols].to_csv(out / "leiden_subset_assignments.csv")
    print(f"\n✓ Leiden assignments saved to {out}/leiden_subset_assignments.csv")

    # ── 5. Marker scoring ────────────────────────────────────────────────────
    print("\nScoring progenitor lineage markers...")
    adata_prog = score_markers(adata_prog)

    # ── 6. UMAPs — QC and batch ──────────────────────────────────────────────
    print("\nSaving UMAPs...")
    for col in [CLUSTER_KEY, "manual_level1", "manual_level2"]:
        save_umap(adata_prog, col, out, "prog_orig_umap")
    sc.pl.embedding(adata_prog, "X_umap_original",
                    color=[CLUSTER_KEY, "manual_level2"], show=False, ncols=2)
    plt.savefig(out / "prog_original_umap_overview.png", dpi=150, bbox_inches="tight")
    plt.close()

    for res in RESOLUTIONS:
        key = f"leiden_sub_r{res}"
        save_umap(adata_prog, key, out, f"r{res}")

    # QC overlays on re-embedded UMAP
    for col in ["pct_counts_mt", "n_genes_by_counts", "n_counts"]:
        if col in adata_prog.obs.columns:
            save_umap(adata_prog, col, out, "qc")
    for col in ["score_Stress_IEG", "score_Cycling"]:
        if col in adata_prog.obs.columns:
            save_umap(adata_prog, col, out, "state")

    # Lineage marker UMAPs
    for lineage, genes in PROGENITOR_MARKERS.items():
        score_col = f"score_{lineage}"
        if score_col in adata_prog.obs.columns:
            save_umap(adata_prog, score_col, out, f"score_{lineage}")
        for g in present(adata_prog, genes[:2]):
            save_umap(adata_prog, g, out, f"gene_{lineage}")

    # ── 7. DE analysis + marker hit tables ───────────────────────────────────
    print("\nRunning DE analysis at each resolution...")
    marker_hits: dict[str, pd.DataFrame] = {}
    for res in RESOLUTIONS:
        key = f"leiden_sub_r{res}"
        print(f"  DE at r={res} ({adata_prog.obs[key].nunique()} clusters vs background)...")
        de_df = de_vs_background(adata_prog, key, out, tag=f"sub_r{res}_vs_background")
        hits = markers_in_top30(de_df)
        hits.to_csv(out / f"marker_hits_r{res}.csv", index=False)
        marker_hits[key] = hits

        # Dotplot per resolution
        save_dotplot(adata_prog, key, out, tag=f"r{res}")

    # DE vs FULL DATASET for each original cluster (critical for IEG-masked clusters)
    print("\nRunning DE for each original progenitor cluster vs full dataset...")
    clean_full_mask = ~pd.Series(adata.var_names).apply(is_problem_gene).values
    a_full_clean = adata[:, clean_full_mask].copy()
    if "lognorm" in a_full_clean.layers:
        a_full_clean.X = a_full_clean.layers["lognorm"].copy()
    sc.tl.rank_genes_groups(
        a_full_clean, groupby=CLUSTER_KEY, method="wilcoxon",
        groups=PROG_CLUSTERS, use_raw=False, key_added="rank_genes_full",
    )
    rows_full = []
    for c in PROG_CLUSTERS:
        names  = a_full_clean.uns["rank_genes_full"]["names"][c][:50]
        scores = a_full_clean.uns["rank_genes_full"]["scores"][c][:50]
        pvals  = a_full_clean.uns["rank_genes_full"]["pvals_adj"][c][:50]
        for rank, (g, s, p) in enumerate(zip(names, scores, pvals)):
            rows_full.append({"cluster_r1.0": c, "rank": rank + 1, "gene": g,
                               "score": round(float(s), 3), "pval_adj": float(p)})
    de_full_df = pd.DataFrame(rows_full)
    de_full_df.to_csv(out / "de_original_clusters_vs_full_dataset.csv", index=False)
    print(f"  ✓ Saved de_original_clusters_vs_full_dataset.csv")
    print("\nTop 20 DE genes per original progenitor cluster (vs full dataset):")
    for c in PROG_CLUSTERS:
        top = de_full_df[de_full_df["cluster_r1.0"] == c].head(20)["gene"].tolist()
        print(f"  Cluster {c}: {top}")

    # ── 8. Annotation template ───────────────────────────────────────────────
    print("\nBuilding annotation template TSV...")
    build_annotation_template(adata_prog, RESOLUTIONS, marker_hits, out)

    # ── 9. Summary table ─────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("PROGENITOR SUBSET — CLUSTER SUMMARY")
    print("=" * 60)
    for res in RESOLUTIONS:
        key = f"leiden_sub_r{res}"
        print(f"\n  Resolution {res}:")
        summary = adata_prog.obs.groupby(key, observed=True).apply(
            lambda g: pd.Series({
                "n_cells": len(g),
                "from_r1.0_clusters": g[CLUSTER_KEY].value_counts().to_dict(),
            })
        )
        print(summary.to_string())
    print("\n✓ All outputs saved to:", out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="results/12_subset_recluster/adata_annotated_final.h5ad",
        help="Path to annotated AnnData",
    )
    parser.add_argument(
        "--out",
        default="results/13_progenitor_recluster",
        help="Output directory",
    )
    main(parser.parse_args())
