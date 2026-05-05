#!/usr/bin/env python
"""Generate final diagnostic figures for the manual cell type annotation."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 1

CORE_MARKERS = [
    "Elane", "Mpo", "Prtn3", "Ms4a3",
    "Mki67", "Top2a",
    "Ltf", "Ngp", "Camp", "Lcn2",
    "Cxcr2", "Ccl6", "Retnlg", "Mmp8",
    "Hlf", "Mecom", "Procr", "Kit", "Ly6a", "Cd34",
    "Ccr2", "Ly6c2", "Csf1r",
    "Cd79a", "Pax5", "Ighm", "Ebf1",
    "Cd3d", "Trac", "Skap1",
    "Irf8", "Siglech", "Bst2", "Cst3",
    "Gata1", "Car2", "Klf1",
    "Pf4", "Itga2b",
    "Ms4a2", "Cpa3", "Gata2",
    "Jchain", "Xbp1", "Mzb1",
]

NEUTROPHIL_MATURATION = [
    "Elane", "Mpo", "Prtn3", "Ms4a3", "Ctsg",
    "Mki67", "Top2a",
    "Ltf", "Ngp", "Camp", "Lcn2", "Cebpe", "Fcnb",
    "Cxcr2", "Ccl6", "Retnlg", "Mmp8", "Mmp9", "S100a8", "S100a9",
]


def safe_genes(adata: sc.AnnData, genes: list[str]) -> list[str]:
    ok = [g for g in genes if g in adata.var_names]
    missing = set(genes) - set(ok)
    if missing:
        print(f"  [warn] not in var_names: {sorted(missing)}")
    return ok


def umap(adata: sc.AnnData, col: str, out: Path, **kw) -> None:
    sc.pl.umap(adata, color=col, show=False, frameon=False,
               legend_loc="right margin", legend_fontsize=7, **kw)
    plt.savefig(out / f"umap_{col}.png", dpi=150, bbox_inches="tight")
    plt.close()


def dotplot(adata: sc.AnnData, genes: list[str], groupby: str, out: Path, fname: str, title: str = "") -> None:
    g = safe_genes(adata, genes)
    if not g:
        return
    sc.pl.dotplot(adata, var_names=g, groupby=groupby, show=False,
                  standard_scale="var", title=title)
    plt.savefig(out / fname, dpi=150, bbox_inches="tight")
    plt.close()


def write_tables(adata: sc.AnnData, out: Path) -> None:
    obs = adata.obs.copy()
    pairs = [
        ("manual_level1", "population", "composition_level1_by_gate.csv"),
        ("manual_level1", "sample_id", "composition_level1_by_sample.csv"),
        ("manual_level2", "population", "composition_level2_by_gate.csv"),
        ("manual_level1", "manual_level2", "level1_vs_level2_crosstab.csv"),
        ("manual_level1", "manual_annotation_confidence", "confidence_by_level1.csv"),
    ]
    for c1, c2, fn in pairs:
        if c1 in obs and c2 in obs:
            pd.crosstab(obs[c1].astype(str), obs[c2].astype(str)).to_csv(out / fn)

    counts = (
        obs.groupby(["manual_level1", "manual_level2"], observed=True)
        .size().reset_index(name="n_cells")
    )
    counts["pct"] = (counts["n_cells"] / adata.n_obs * 100).round(2)
    counts.to_csv(out / "level1_level2_cell_counts.csv", index=False)


def write_attestation(adata: sc.AnnData, out: Path) -> None:
    text = f"""\
# Annotation Attestation — Charles Mouse BM scRNA-seq
# leiden_r1.0, {adata.n_obs} cells, {len(adata.obs['manual_level1'].unique())} level1 categories

## Resolution
Annotation performed at leiden_r1.0 (23 clusters).
leiden_r2.0 was evaluated and rejected: 3 clusters < 100 cells (min=23 cells),
6 clusters < 300 cells. The two meaningful splits at r2.0 (B cell pro/pre-B vs
mature-B; erythroid MEP vs erythroblast) do not justify the noise floor.

## What Was Annotated with Confidence

### High confidence, both levels
- granulocyte_neutrophil:
    - GMP_neutrophil_primed (clusters 2, 5): LK gate, Elane/Mpo/Prtn3/Ms4a3
    - Immature_neutrophil (clusters 4, 9, 14, 16): Ltf/Ngp/Camp/Lcn2
    - Mature_neutrophil (clusters 0, 1, 11): Cxcr2/Ccl6/Mmp8/Retnlg; I gate
- cycling_myeloid (clusters 7, 8): Mki67/Top2a + LK gate
- monocyte / Classical_monocyte (cluster 12): Ccr2/Ly6c2/Csf1r
- DC / pDC (cluster 20): Siglech/Bst2/Tcf4 — highest-confidence cluster in dataset
- plasma_cell (cluster 22): Jchain/Mzb1/Xbp1

### Medium confidence
- DC / DC_broad (cluster 17): Irf8/Ifi205/Id2 — cDC1-biased; Clec9a/Xcr1 absent
- erythroid / Erythroid_broad (cluster 13): Car2/Klf1 erythroid program clear;
  MEP vs Erythroblast substage NOT resolvable at this resolution
- B_cell / B_cell_broad (cluster 15): Cd79a/Ighm/Pax5 clear; Pro-Pre-B vs Mature-B
  NOT resolvable (Vpreb3 and H2-Aa both present in single cluster)
- granulocyte_neutrophil / Mature_neutrophil (cluster 11): Il1b/Malat1 co-expression
  may indicate inflammatory or low-quality subpopulation; flag for sensitivity analysis
- granulocyte_neutrophil / Mature_neutrophil (cluster 21): 78 cells, low statistical weight

## What CANNOT Be Resolved — Do Not Report These in the Manuscript

1. T cell subtype (CD4/CD8/γδ): 1 cluster, 314 cells. CD4/CD8 not in top DE genes.
   ANNOTATION: T_cell / T_cell_broad (level2 = T_cell_broad is a placeholder, not a finding)

2. NK cells: not separated from T cells at this clustering resolution.
   Nkg7/Klrb1c partially present in cluster 18 but insufficient for confident NK label.

3. B cell pro/pre-B vs mature-B: Vpreb3 (pro-B marker) and H2-Aa (mature B marker)
   co-present in cluster 15 (807 cells). Resolution requires subset reclustering or r2.0.

4. Erythroid MEP vs Erythroblast: Car2/Blvrb and early markers co-present in cluster 13.
   Resolution requires subset reclustering or r2.0.

5. HSC vs MPP vs LMPP: clusters 3 and 10 both show high mitochondrial gene content
   (mt-Co1/mt-Co3/mt-Atp8 in top DE), confounding marker scoring.
   Cannot separate HSC subtypes without RNA velocity or surface protein co-staining.

6. Cluster 6 (myeloid_progenitor / Progenitor_ambiguous): popV assigns CMP(59%)/MEP(22%)/
   Basophil(6%). Gata2 is present but basophil level1 score = 0.66 (insufficient).
   Basophil progenitors are likely a subpopulation but cannot be isolated at r1.0.

7. DC cDC1 vs cDC2: cluster 17 (378 cells) shows cDC1-biased markers (Ifi205/Id2)
   but Clec9a/Xcr1 are absent from top DE genes. Classified DC_broad.

8. Cluster 19 (119 cells): Dock2/Elmo1/Arhgap15 lymphoid-migration program in I gate.
   No canonical lineage can be assigned.

## Recommended Next Steps for Higher Granularity
- Targeted subset reclustering of cluster 6 (myeloid_progenitor) to isolate basophil progenitors
- Targeted subset reclustering of cluster 15 (B_cell) to separate pro/pre-B from mature-B
- Targeted subset reclustering of cluster 13 (erythroid) to separate MEP from erythroblast
- RNA velocity (scVelo) on progenitor + neutrophil compartment for pseudotime ordering
- CITE-seq or FlowJo surface staining for CD4/CD8 to subtype T cells
"""
    (out / "ANNOTATION_ATTESTATION.md").write_text(text)


def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"]

    for col in ["manual_level1", "manual_level2", "manual_annotation_confidence",
                "manual_review_flag", "leiden_r1.0", "population"]:
        if col in adata.obs:
            print(f"  UMAP: {col}")
            umap(adata, col, out)

    dotplot(adata, CORE_MARKERS, "manual_level1", out,
            "dotplot_core_markers_by_level1.png", "Core markers by level1")

    dotplot(adata, NEUTROPHIL_MATURATION, "leiden_r1.0", out,
            "dotplot_neutrophil_maturation_by_cluster.png",
            "Neutrophil maturation gradient by cluster (r1.0)")

    neu = adata[adata.obs["manual_level1"].isin(
        ["granulocyte_neutrophil", "myeloid_progenitor", "cycling_myeloid"]
    )].copy()
    if neu.n_obs > 50:
        dotplot(neu, NEUTROPHIL_MATURATION, "manual_level2", out,
                "dotplot_neutrophil_maturation_by_level2.png",
                "Neutrophil maturation by level2 (myeloid compartment)")

    write_tables(adata, out)
    write_attestation(adata, out)
    print(f"Done. Outputs in {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/12_subset_recluster/adata_annotated_final.h5ad")
    parser.add_argument("--out", default="results/13_final_figures")
    main(parser.parse_args())
