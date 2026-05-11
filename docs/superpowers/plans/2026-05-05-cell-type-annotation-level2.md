# Cell Type Annotation — Level2 Completion Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete coarse-to-fine cell type annotation of the Charles mouse bone marrow scRNA-seq dataset by adding level2 labels to all 23 `leiden_r1.0` clusters, using sorted gate information as the primary biological constraint and DE gene evidence for level2 calls.

**Architecture:** Level1 annotation was previously applied at `leiden_r2.0` resolution (`results/10_manual_level1/`). This plan switches to `leiden_r1.0` (23 clusters) as the primary annotation resolution, creates a complete gate-constrained level1+level2 map as a static TSV, re-runs script 10 with the r1.0 cluster key, then performs targeted subset reclustering for the two clusters where level2 is not resolvable at global r1.0 resolution (B cell cluster 15 and erythroid cluster 13), and finally generates final diagnostic figures.

**Tech Stack:** Python 3.11 + Scanpy + Matplotlib via `pixi run python`, TSV config files.

---

## Resolution Decision: Why leiden_r1.0 and Not r2.0

`leiden_r2.0` was evaluated and rejected as the primary annotation resolution.

**Evidence against r2.0:**

| Issue | Data |
|-------|------|
| Cluster 35: 23 cells | Too few for reliable DE and marker scoring |
| Cluster 33: 74 cells, cluster 31: 199 cells | Marginally informative |
| 6 of 36 clusters < 300 cells | Exceeds acceptable noise floor for this data quality |
| r2.0 neutrophil fragmentation | 12 neutrophil clusters at r1.0 → 20+ at r2.0, same DE profile |
| Progenitor over-splitting | r1.0 clusters 3+10 split into 4 r2.0 clusters, only 1 (cluster 17) gains identity |

**What r2.0 adds over r1.0 (the genuine splits):**

| Population | r1.0 | r2.0 | Gain |
|------------|------|------|------|
| B cell | cluster 15, 807 cells | clusters 26+27 (Vpreb3 vs H2-Aa) | Pro-Pre-B vs Mature-B |
| Erythroid | cluster 13, 1062 cells | clusters 13+24 (MEP vs Erythroblast) | MEP vs Erythroblast |

These two splits are biologically meaningful but do not justify accepting 3 clusters with < 100 cells and 6 with < 300 cells across the full object. The appropriate strategy for these two compartments is targeted subset analysis (optional extension, not required for the first annotation pass).

**Conclusion: use `leiden_r1.0` (23 clusters) as the annotation basis.**

---

## Background: Sorted Gate Constraints

The sorted fractions serve as hard biological priors. No annotation can contradict them.

| Gate | Biological identity | Dominant clusters (r1.0) |
|------|---------------------|--------------------------|
| `LSK` (Lin− Sca1+ Kit+) | HSC, MPP, LMPP | 3, 10 (mostly) |
| `LK` (Lin− Kit+ Sca1−) | GMP, CMP, MEP, immature granulocyte | 2, 5, 6, 7, 8, 13, 14, 17 |
| `I` (Input, unsorted) | Mature neutrophil, monocyte, B/T cell, pDC, plasma | 0, 1, 4, 9, 11, 12, 15, 16, 18, 19, 20, 21, 22 |

---

## Annotation Ceiling: What Can and Cannot Be Resolved

**Achievable at r1.0:**
- Neutrophil maturation gradient: GMP/Pro-neutrophil → Immature → Mature
- Cycling myeloid (Mki67+/Top2a+ in LK gate)
- Classical monocyte (Ccr2/Ly6c2, I gate)
- DC vs pDC (already separated at r1.0: clusters 17 and 20)
- B cell (single cluster — broad call only; see below)
- T cell (single cluster)
- Erythroid (single cluster — broad call only)
- Megakaryocyte, plasma cell: embedded in mixed cluster 6 and cluster 22 respectively

**NOT achievable at this depth — formally attested in output:**
- T cell subtype (CD4/CD8/γδ): 1 cluster, 314 cells, CD4/CD8 not in top DE genes
- NK cells: not separated from T cells; Nkg7/Klrb1c partially overlap with cluster 18
- B cell pro/pre-B vs mature-B: r1.0 cluster 15 mixes both (Vpreb3 and H2-Aa co-present); requires targeted subset or r2.0 upgrade for this compartment
- Erythroid MEP vs Erythroblast: r1.0 cluster 13 mixes both; Car2/Klf1 and early markers co-present
- HSC vs MPP vs LMPP: clusters 3 and 10 both have high mitochondrial content in top DE, confounding marker scoring
- Basophil vs CMP/MEP within cluster 6: heterogeneous cluster (popV: CMP 59%, MEP 22%, basophil 6%); basophil level1 score = 0.66 (too weak for confident level1 basophil call)
- DC cDC1 vs cDC2: cluster 17 shows Ifi205/Id2 (cDC1-biased) but Clec9a/Xcr1 absent from top DE
- Ambiguous cluster 19 (119 cells): lymphoid-migration DE (Dock2/Elmo1/Arhgap15) in I gate; no canonical lineage

---

## Complete Level2 Annotation Decision Table (leiden_r1.0, 23 clusters)

Decision logic for the neutrophil compartment:
- `GMP_neutrophil_primed`: LK gate dominant + Elane/Mpo/Prtn3/Ctsg high + Ltf/Ngp/Cxcr2 absent
- `Immature_neutrophil`: Ltf+/Ngp+/Camp+/Lcn2+ predominant, Cxcr2 absent (any gate)
- `Mature_neutrophil`: Cxcr2+/Ccl6+/Mmp8+/Retnlg+ predominant, I gate dominant
- `Cycling_myeloid`: Mki67+/Top2a+ in LK gate context

| cluster | n_cells | gate | level1 | level2 | confidence | key DE evidence |
|---------|---------|------|--------|--------|------------|-----------------|
| 0 | 6084 | I | granulocyte_neutrophil | Mature_neutrophil | high | Mmp8, Retnlg, Mmp9, S100a8, S100a9, S100a6 |
| 1 | 4791 | I | granulocyte_neutrophil | Mature_neutrophil | high | Ccl6, Tyrobp, Csf3r, Cxcr2, Clec4d |
| 2 | 4291 | LK | myeloid_progenitor | GMP_neutrophil_primed | high | Prtn3, Mpo, Elane, Ctsg, Plac8 — azurophil granule program, no Ltf/Cxcr2 |
| 3 | 3963 | LSK | progenitor | Progenitor_ambiguous | low | Adgrl4, Cd34, Angpt1, mt-Co1 — HSPC markers present but high mito contaminates |
| 4 | 3749 | I | granulocyte_neutrophil | Immature_neutrophil | high | Ltf, Anxa1, Ngp, Lcn2, Camp, Chil3 |
| 5 | 2803 | LK | myeloid_progenitor | GMP_neutrophil_primed | high | Mpo, Ctsg, Elane, Prtn3, Ms4a3 — Ms4a3 = high-confidence GMP marker |
| 6 | 2691 | LK | myeloid_progenitor | Progenitor_ambiguous | low | Apoe, Gata2, Jun, Angpt1, Egr1, Cdk6 — heterogeneous: popV=CMP(59%)/MEP(22%)/Bas(6%); basophil score 0.66; cannot resolve without subset reclustering |
| 7 | 2641 | LK | cycling_myeloid | Cycling_myeloid | high | Tubb5, Tuba1b, Birc5, H2afx, Pclaf, Aurkb, Top2a |
| 8 | 1670 | LK | cycling_myeloid | Cycling_myeloid | high | Top2a, Mki67, Smc4, Kif11, Knl1 |
| 9 | 1542 | I | granulocyte_neutrophil | Immature_neutrophil | high | S100a9, Camp, S100a8, Ngp, Wfdc21, Ifitm6, Lcn2 |
| 10 | 1297 | LSK | progenitor | Progenitor_ambiguous | low | Angpt1, mt-Co1, mt-Co3, mt-Co2, mt-Atp8 — dominated by mitochondrial genes |
| 11 | 1158 | I | granulocyte_neutrophil | Mature_neutrophil | medium | Csf3r, Il1b, Tyrobp, Malat1 — Il1b/Malat1 suggest inflammatory/low-quality state |
| 12 | 1096 | I | monocyte | Classical_monocyte | high | Ctss, S100a4, Psap, Ccr2, Ahnak, Pld4 — Ccr2/Ly6c2 canonical classical monocyte |
| 13 | 1062 | LK | erythroid | Erythroid_broad | medium | Car1, Blvrb, Car2, Klf1, Atpif1 — clean erythroid program; MEP vs Erythroblast split NOT resolvable at r1.0 |
| 14 | 1041 | I | granulocyte_neutrophil | Immature_neutrophil | high | Hmgn2, Chil3, Camp, Fcnb, Hmgb2, Cebpe — Fcnb=neutrophil granule protein |
| 15 | 807 | I | B_cell | B_cell_broad | medium | Cd79a, Ighm, Ebf1, Bach2, Vpreb3, Pax5 — Pro-Pre-B vs Mature-B NOT resolvable at r1.0 (Vpreb3 and H2-Aa both present in the cluster) |
| 16 | 638 | I | granulocyte_neutrophil | Immature_neutrophil | high | Chil3, Camp, Ngp, S100a8, Orm1, Lcn2 |
| 17 | 378 | LK | DC | DC_broad | medium | Cst3, Irf8, Ifi205, Cd74, Id2, H2-Aa — Ifi205/Id2 suggest cDC1-biased; Clec9a/Xcr1 absent from top DE |
| 18 | 314 | I | T_cell | T_cell_broad | high | Tmsb10, Skap1, Ms4a4b, H2-Q7 — T cell clear; CD4/CD8 NOT in top DE; subtype unresolvable |
| 19 | 119 | I | ambiguous_myeloid | Ambiguous | low | Malat1, Arhgap15, Dock2, Baz2b, Elmo1 — lymphoid-migration program in I gate; lineage unresolvable |
| 20 | 115 | I | DC | pDC | high | Tcf4, Irf8, Bst2, Siglech — definitive pDC markers |
| 21 | 78 | I | granulocyte_neutrophil | Mature_neutrophil | low | Lgals3, Pfn1, Mmp8, Msrb1 — 78 cells; Mmp8 + I gate = mature; small cluster, low weight |
| 22 | 62 | I/LSK | plasma_cell | Plasma_cell | high | Jchain, Mzb1, Xbp1, Igkc, Txndc5 — all canonical plasma cell markers |

---

## Task 1: Write the complete annotation map TSV for leiden_r1.0

**Files:**
- Create: `pipeline/config/manual_annotation_full_map_leiden_r1.0.tsv`

- [ ] **Step 1: Create the TSV**

```bash
cat > pipeline/config/manual_annotation_full_map_leiden_r1.0.tsv << 'EOF'
cluster_key	cluster	manual_level1	manual_level2	manual_annotation_confidence	manual_review_flag	manual_annotation_basis
leiden_r1.0	0	granulocyte_neutrophil	Mature_neutrophil	high	ok	Mmp8/Retnlg/Mmp9/S100a8/S100a9; I gate dominant
leiden_r1.0	1	granulocyte_neutrophil	Mature_neutrophil	high	ok	Ccl6/Tyrobp/Csf3r/Cxcr2/Clec4d; I gate dominant
leiden_r1.0	2	myeloid_progenitor	GMP_neutrophil_primed	high	ok	Prtn3/Mpo/Elane/Ctsg/Plac8; azurophil granule program; Ltf/Ngp/Cxcr2 absent; LK gate
leiden_r1.0	3	progenitor	Progenitor_ambiguous	low	high_mito_review	Adgrl4/Cd34/Angpt1 HSPC-like but mt-Co1 in top DE; LSK gate dominant; cannot subtype HSC/MPP/LMPP
leiden_r1.0	4	granulocyte_neutrophil	Immature_neutrophil	high	ok	Ltf/Anxa1/Ngp/Lcn2/Camp/Chil3; I gate dominant
leiden_r1.0	5	myeloid_progenitor	GMP_neutrophil_primed	high	ok	Mpo/Ctsg/Elane/Prtn3/Ms4a3; Ms4a3=high-confidence GMP marker; LK gate
leiden_r1.0	6	myeloid_progenitor	Progenitor_ambiguous	low	ambiguous_cmp_mep_basophil	Gata2/Angpt1/Apoe; popV=CMP(59%)/MEP(22%)/Bas(6%); basophil level1 score 0.66 too weak; heterogeneous LK cluster; cannot resolve without subset reclustering
leiden_r1.0	7	cycling_myeloid	Cycling_myeloid	high	ok	Tubb5/Tuba1b/Birc5/H2afx/Pclaf/Aurkb/Top2a; LK gate
leiden_r1.0	8	cycling_myeloid	Cycling_myeloid	high	ok	Top2a/Mki67/Smc4/Kif11/Knl1; LK gate
leiden_r1.0	9	granulocyte_neutrophil	Immature_neutrophil	high	ok	S100a9/Camp/S100a8/Ngp/Wfdc21/Ifitm6/Lcn2; I gate dominant
leiden_r1.0	10	progenitor	Progenitor_ambiguous	low	high_mito_review	mt-Co1/mt-Co3/mt-Co2/mt-Atp8 dominate top DE; mitochondrial content too high for progenitor subtyping; LSK gate
leiden_r1.0	11	granulocyte_neutrophil	Mature_neutrophil	medium	inflammatory_review	Csf3r/Il1b/Tyrobp/Malat1; Il1b and Malat1 suggest inflammatory or low-quality mature neutrophil; I gate
leiden_r1.0	12	monocyte	Classical_monocyte	high	ok	Ctss/S100a4/Psap/Ccr2/Ahnak/Pld4; Ccr2/Ly6c2 canonical classical monocyte; I gate
leiden_r1.0	13	erythroid	Erythroid_broad	medium	mep_vs_erythroblast_unresolvable	Car1/Blvrb/Car2/Klf1/Atpif1; clean erythroid program but MEP vs Erythroblast split not resolvable at r1.0 granularity (single cluster)
leiden_r1.0	14	granulocyte_neutrophil	Immature_neutrophil	high	ok	Hmgn2/Chil3/Camp/Fcnb/Hmgb2/Cebpe; Fcnb=neutrophil granule; I gate dominant
leiden_r1.0	15	B_cell	B_cell_broad	medium	b_cell_subtype_unresolvable	Cd79a/Ighm/Ebf1/Bach2/Vpreb3/Pax5; Vpreb3 (pro-B) and H2-Aa (mature-B) co-present; Pro-Pre-B vs Mature-B NOT resolvable at r1.0 (single cluster, 807 cells)
leiden_r1.0	16	granulocyte_neutrophil	Immature_neutrophil	high	ok	Chil3/Camp/Ngp/S100a8/Orm1/Lcn2; I gate dominant
leiden_r1.0	17	DC	DC_broad	medium	cdc1_cdc2_review	Cst3/Irf8/Ifi205/Cd74/Id2/H2-Aa; Ifi205/Id2 suggest cDC1-biased but Clec9a/Xcr1 absent from top DE; classify DC_broad
leiden_r1.0	18	T_cell	T_cell_broad	high	subtype_unresolvable	Tmsb10/Skap1/Ms4a4b/H2-Q7; T cell clear but CD4/CD8 absent from top DE; T cell subtype CANNOT be resolved at this depth
leiden_r1.0	19	ambiguous_myeloid	Ambiguous	low	lineage_unresolvable	Malat1/Arhgap15/Dock2/Baz2b/Elmo1; lymphoid-migration program in I gate; cannot assign canonical lineage
leiden_r1.0	20	DC	pDC	high	ok	Tcf4/Irf8/Bst2/Siglech; definitive pDC markers; I gate
leiden_r1.0	21	granulocyte_neutrophil	Mature_neutrophil	low	small_cluster	Lgals3/Pfn1/Mmp8/Msrb1; 78 cells; Mmp8+I gate=mature; small cluster with low interpretive weight
leiden_r1.0	22	plasma_cell	Plasma_cell	high	ok	Jchain/Mzb1/Xbp1/Igkc/Txndc5; canonical plasma cell; I gate
EOF
```

- [ ] **Step 2: Validate column count (must be 7 for all rows)**

```bash
awk -F '\t' 'NF!=7 {print NR, NF, $0}' pipeline/config/manual_annotation_full_map_leiden_r1.0.tsv
```

Expected: no output.

- [ ] **Step 3: Verify all 23 clusters present**

```bash
tail -n +2 pipeline/config/manual_annotation_full_map_leiden_r1.0.tsv | cut -f2 | sort -n | tr '\n' ','
```

Expected: `0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,`

- [ ] **Step 4: Commit the map**

```bash
git add pipeline/config/manual_annotation_full_map_leiden_r1.0.tsv
git commit -m "annotation: add complete level1+level2 map for leiden_r1.0 (gate-constrained, DE-guided)"
```

---

## Task 2: Patch script 10 to apply level2 from the map

**Files:**
- Modify: `pipeline/scripts/10_apply_manual_level1_annotation.py`

Script 10 currently hardcodes `manual_level2 = "not_assigned_level1_only"` (line ~76). The map now has a `manual_level2` column; the script must read it.

- [ ] **Step 1: Confirm the hardcoded line**

```bash
grep -n "not_assigned_level1_only" pipeline/scripts/10_apply_manual_level1_annotation.py
```

Expected: one match around line 76.

- [ ] **Step 2: Apply Edit A — add `manual_level2` to required columns**

In `read_annotation_map`, around line 24, update `required`:

```python
# Before:
required = {
    "cluster_key",
    "cluster",
    "manual_level1",
    "manual_annotation_confidence",
    "manual_review_flag",
    "manual_annotation_basis",
}

# After:
required = {
    "cluster_key",
    "cluster",
    "manual_level1",
    "manual_level2",
    "manual_annotation_confidence",
    "manual_review_flag",
    "manual_annotation_basis",
}
```

- [ ] **Step 3: Apply Edit B — replace hardcoded level2 with map lookup in `apply_map`**

```python
# Before (around line 76):
adata.obs["manual_level2"] = "not_assigned_level1_only"

# After:
adata.obs["manual_level2"] = clusters.map(keyed["manual_level2"]).astype("category").to_numpy()
```

- [ ] **Step 4: Apply Edit C — add `manual_level2` to cluster summary rows in `write_cluster_summary`**

In the `row` dict inside `write_cluster_summary` (around line 95), add after `manual_level1`:

```python
"manual_level2": sub["manual_level2"].astype(str).iloc[0],
```

- [ ] **Step 5: Apply Edit D — add `manual_level2` UMAP to `plot_umaps`**

```python
# Before:
for col in ["manual_level1", "manual_annotation_confidence", "manual_review_flag"]:

# After:
for col in ["manual_level1", "manual_level2", "manual_annotation_confidence", "manual_review_flag"]:
```

- [ ] **Step 6: Verify the patch compiles**

```bash
pixi run python -m py_compile pipeline/scripts/10_apply_manual_level1_annotation.py && echo "OK"
```

Expected: `OK`

- [ ] **Step 7: Commit**

```bash
git add pipeline/scripts/10_apply_manual_level1_annotation.py
git commit -m "feat: patch script 10 to read manual_level2 from annotation map TSV"
```

---

## Task 3: Run annotation at leiden_r1.0 with the full map

**Files:**
- Input: `results/08_reconcile_annotations/adata_reconciled.h5ad`
- Config: `pipeline/config/manual_annotation_full_map_leiden_r1.0.tsv`
- Output: `results/11_annotation_r1.0/`

- [ ] **Step 1: Run annotation**

```bash
pixi run python pipeline/scripts/10_apply_manual_level1_annotation.py \
    --input results/08_reconcile_annotations/adata_reconciled.h5ad \
    --map pipeline/config/manual_annotation_full_map_leiden_r1.0.tsv \
    --cluster-key leiden_r1.0 \
    --annotation-version manual_level1_level2_r1.0_2026-05-05 \
    --out results/11_annotation_r1.0
```

Expected: `Saved manual level1 annotation outputs to results/11_annotation_r1.0`

- [ ] **Step 2: Verify level2 is populated and not the hardcoded placeholder**

```bash
pixi run python -c "
import pandas as pd
df = pd.read_csv('results/11_annotation_r1.0/manual_level1_per_cell.csv', index_col=0)
print('level2 unique values:')
print(df['manual_level2'].value_counts().to_string())
print()
print('Hardcoded placeholder present:', df['manual_level2'].eq('not_assigned_level1_only').any())
"
```

Expected: `Hardcoded placeholder present: False` and a distribution showing Mature_neutrophil, Immature_neutrophil, GMP_neutrophil_primed, etc.

- [ ] **Step 3: Verify cluster summary looks correct**

```bash
pixi run python -c "
import pandas as pd
df = pd.read_csv('results/11_annotation_r1.0/manual_level1_cluster_summary_leiden_r1.0.csv')
print(df[['cluster','manual_level1','manual_level2','n_cells','manual_annotation_confidence']].to_string(index=False))
"
```

Expected: 23 rows, no `not_assigned_level1_only` in level2 column.

- [ ] **Step 4: Commit tabular outputs**

```bash
git add results/11_annotation_r1.0/*.csv results/11_annotation_r1.0/*.tsv
git commit -m "results: level1+level2 annotation applied at leiden_r1.0 (23 clusters)"
```

---

## Task 4: Write subset reclustering script for B cell and erythroid compartments

**Files:**
- Create: `pipeline/scripts/12_subset_recluster.py`

At global `leiden_r1.0`, two clusters cannot be resolved at level2:
- Cluster 15 (807 cells, B cell): Vpreb3 (pro/pre-B) and H2-Aa (mature-B) co-present in a single cluster
- Cluster 13 (1062 cells, erythroid): MEP markers (Gata1/Klf1) and erythroblast markers (Car2/Blvrb) co-present

This script extracts each cluster from the annotated object, re-embeds with PCA→neighbors→leiden at subset-level, scores lineage markers, assigns level2 labels, and updates `manual_level2` in the full object. If the subset reclustering does not produce a meaningful split (i.e., both marker scores point to the same call), the script keeps the broad label unchanged and reports this in a summary CSV.

Input: `results/11_annotation_r1.0/adata_manual_level1.h5ad`
Output: `results/12_subset_recluster/adata_annotated_final.h5ad` + diagnostic plots + summary CSVs

- [ ] **Step 1: Create the script**

```python
#!/usr/bin/env python
"""Subset reclustering for B cell and erythroid compartments.

For each compartment, extract the relevant r1.0 cluster, re-embed at higher
resolution, score lineage markers, assign level2 labels, and update the full
object's manual_level2 column. If the split is not supported by the data,
keep the broad label and record the outcome in a CSV.
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

B_CLUSTER = "15"
ERY_CLUSTER = "13"
CLUSTER_KEY = "leiden_r1.0"

# Markers: B cell
PRO_PRE_B = ["Vpreb1", "Vpreb3", "Ebf1", "Dntt", "Bach2"]
MATURE_B   = ["H2-Aa", "H2-Eb1", "H2-Ab1", "Bank1", "Ms4a1", "Cd74"]

# Markers: erythroid
MEP_ERY    = ["Gata1", "Klf1", "Itga2b", "Gata2", "Hba-a1"]
ERYTHROBLAST = ["Car2", "Car1", "Blvrb", "Aqp1", "Ermap", "Slc4a1"]


def present(adata: sc.AnnData, genes: list[str]) -> list[str]:
    return [g for g in genes if g in adata.var_names]


def recluster(adata_sub: sc.AnnData, res: float = 0.5, n_hvg: int = 400,
              n_pcs: int = 15, n_neighbors: int = 10) -> sc.AnnData:
    """Re-embed a subset: HVG → PCA → neighbors → leiden → UMAP."""
    if "lognorm" in adata_sub.layers:
        adata_sub.X = adata_sub.layers["lognorm"].copy()
    sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_hvg, subset=False)
    sc.pp.pca(adata_sub, n_comps=n_pcs, use_highly_variable=True)
    sc.pp.neighbors(adata_sub, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(adata_sub, resolution=res, key_added="leiden_subset")
    sc.tl.umap(adata_sub)
    return adata_sub


def mean_expr(adata_sub: sc.AnnData, genes: list[str]) -> float:
    ok = present(adata_sub, genes)
    if not ok:
        return 0.0
    return float(np.asarray(adata_sub[:, ok].X.mean()))


def subcluster_calls(adata_sub: sc.AnnData,
                     pos_markers: list[str], neg_markers: list[str],
                     pos_label: str, neg_label: str) -> dict[str, str]:
    """Return {subcluster_id: level2_label} for each leiden_subset cluster."""
    calls: dict[str, str] = {}
    for c in sorted(adata_sub.obs["leiden_subset"].unique(), key=int):
        sub = adata_sub[adata_sub.obs["leiden_subset"] == c]
        pos_score = mean_expr(sub, pos_markers)
        neg_score = mean_expr(sub, neg_markers)
        calls[c] = pos_label if pos_score >= neg_score else neg_label
    return calls


def plot_subset(adata_sub: sc.AnnData, color_cols: list[str],
                out: Path, prefix: str) -> None:
    for col in color_cols:
        if col not in adata_sub.obs and col not in adata_sub.var_names:
            continue
        sc.pl.umap(adata_sub, color=col, show=False, frameon=False,
                   legend_loc="right margin", legend_fontsize=7)
        plt.savefig(out / f"{prefix}_umap_{col.replace('/', '_')}.png",
                    dpi=150, bbox_inches="tight")
        plt.close()


def process_b_cell(adata: sc.AnnData, out: Path) -> pd.DataFrame:
    """Recluster B cell cluster and assign Pro_Pre_B vs Mature_B level2."""
    mask = adata.obs[CLUSTER_KEY].astype(str) == B_CLUSTER
    adata_b = adata[mask].copy()
    print(f"  B cell subset: {adata_b.n_obs} cells")

    adata_b = recluster(adata_b, res=0.5)
    calls = subcluster_calls(adata_b, PRO_PRE_B, MATURE_B, "Pro_Pre_B", "Mature_B")

    rows = []
    for c, label in calls.items():
        sub = adata_b[adata_b.obs["leiden_subset"] == c]
        rows.append({
            "compartment": "B_cell",
            "subcluster": c,
            "n_cells": int(sub.n_obs),
            "pro_pre_b_markers_mean": round(mean_expr(sub, PRO_PRE_B), 4),
            "mature_b_markers_mean": round(mean_expr(sub, MATURE_B), 4),
            "level2_call": label,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "b_cell_subcluster_summary.csv", index=False)

    unique = summary["level2_call"].nunique()
    if unique == 1:
        print("  [B cell] No meaningful split — keeping B_cell_broad")
        adata_b.obs["level2_updated"] = "B_cell_broad"
    else:
        adata_b.obs["level2_updated"] = adata_b.obs["leiden_subset"].map(calls).fillna("B_cell_broad")
        print(f"  [B cell] Split: {dict(summary.groupby('level2_call')['n_cells'].sum())}")

    plot_subset(adata_b, ["leiden_subset", "level2_updated", "Vpreb3", "H2-Aa"], out, "b_cell")
    return adata_b.obs[["level2_updated"]].rename(columns={"level2_updated": "manual_level2_new"})


def process_erythroid(adata: sc.AnnData, out: Path) -> pd.DataFrame:
    """Recluster erythroid cluster and assign MEP_Erythroid vs Erythroblast level2."""
    mask = adata.obs[CLUSTER_KEY].astype(str) == ERY_CLUSTER
    adata_e = adata[mask].copy()
    print(f"  Erythroid subset: {adata_e.n_obs} cells")

    adata_e = recluster(adata_e, res=0.5)
    calls = subcluster_calls(adata_e, MEP_ERY, ERYTHROBLAST, "MEP_Erythroid", "Erythroblast")

    rows = []
    for c, label in calls.items():
        sub = adata_e[adata_e.obs["leiden_subset"] == c]
        rows.append({
            "compartment": "erythroid",
            "subcluster": c,
            "n_cells": int(sub.n_obs),
            "mep_markers_mean": round(mean_expr(sub, MEP_ERY), 4),
            "erythroblast_markers_mean": round(mean_expr(sub, ERYTHROBLAST), 4),
            "level2_call": label,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "erythroid_subcluster_summary.csv", index=False)

    unique = summary["level2_call"].nunique()
    if unique == 1:
        print("  [Erythroid] No meaningful split — keeping Erythroid_broad")
        adata_e.obs["level2_updated"] = "Erythroid_broad"
    else:
        adata_e.obs["level2_updated"] = adata_e.obs["leiden_subset"].map(calls).fillna("Erythroid_broad")
        print(f"  [Erythroid] Split: {dict(summary.groupby('level2_call')['n_cells'].sum())}")

    plot_subset(adata_e, ["leiden_subset", "level2_updated", "Gata1", "Car2"], out, "erythroid")
    return adata_e.obs[["level2_updated"]].rename(columns={"level2_updated": "manual_level2_new"})


def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells")

    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()

    # Run subset reclustering for each compartment
    b_updates   = process_b_cell(adata, out)
    ery_updates = process_erythroid(adata, out)

    # Merge updates back into the full object
    all_updates = pd.concat([b_updates, ery_updates])
    overlap = all_updates.index[all_updates.index.duplicated()]
    if len(overlap):
        raise ValueError(f"Duplicate barcodes in subset updates: {overlap[:5].tolist()}")

    adata.obs["manual_level2"] = adata.obs["manual_level2"].astype(str)
    adata.obs.loc[all_updates.index, "manual_level2"] = all_updates["manual_level2_new"].values
    adata.obs["manual_level2"] = adata.obs["manual_level2"].astype("category")

    # Write final counts
    counts = (
        adata.obs.groupby(["manual_level1", "manual_level2"], observed=True)
        .size().reset_index(name="n_cells")
    )
    counts["pct"] = (counts["n_cells"] / adata.n_obs * 100).round(2)
    counts.to_csv(out / "final_level1_level2_counts.csv", index=False)

    adata.write_h5ad(out / "adata_annotated_final.h5ad")
    print(f"Saved final annotated object to {out}/adata_annotated_final.h5ad")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  default="results/11_annotation_r1.0/adata_manual_level1.h5ad")
    parser.add_argument("--out",    default="results/12_subset_recluster")
    main(parser.parse_args())
```

Save to `pipeline/scripts/12_subset_recluster.py`.

- [ ] **Step 2: Verify compilation**

```bash
pixi run python -m py_compile pipeline/scripts/12_subset_recluster.py && echo "OK"
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add pipeline/scripts/12_subset_recluster.py
git commit -m "feat: subset reclustering script for B cell and erythroid level2 resolution"
```

---

## Task 5: Run subset reclustering and verify level2 updates

**Files:**
- Input: `results/11_annotation_r1.0/adata_manual_level1.h5ad`
- Output: `results/12_subset_recluster/`

- [ ] **Step 1: Run the script**

```bash
pixi run python pipeline/scripts/12_subset_recluster.py \
    --input results/11_annotation_r1.0/adata_manual_level1.h5ad \
    --out results/12_subset_recluster
```

Expected output: lines reporting B cell (807 cells) and erythroid (1062 cells) subset sizes, split outcome for each, and `Saved final annotated object to results/12_subset_recluster/adata_annotated_final.h5ad`

- [ ] **Step 2: Verify split outcomes**

```bash
pixi run python -c "
import pandas as pd
b = pd.read_csv('results/12_subset_recluster/b_cell_subcluster_summary.csv')
e = pd.read_csv('results/12_subset_recluster/erythroid_subcluster_summary.csv')
print('=== B cell subclusters ===')
print(b[['subcluster','n_cells','pro_pre_b_markers_mean','mature_b_markers_mean','level2_call']].to_string(index=False))
print()
print('=== Erythroid subclusters ===')
print(e[['subcluster','n_cells','mep_markers_mean','erythroblast_markers_mean','level2_call']].to_string(index=False))
"
```

- [ ] **Step 3: Verify level2 updates in the final object**

```bash
pixi run python -c "
import pandas as pd
df = pd.read_csv('results/12_subset_recluster/final_level1_level2_counts.csv')
lymph = df[df['manual_level1'].isin(['B_cell','erythroid'])]
print(lymph.to_string(index=False))
"
```

Expected: B_cell rows showing Pro_Pre_B and/or Mature_B (or B_cell_broad if split not supported); erythroid rows showing MEP_Erythroid and/or Erythroblast (or Erythroid_broad if split not supported).

- [ ] **Step 4: Commit outputs**

```bash
git add results/12_subset_recluster/*.csv
git commit -m "results: B cell and erythroid subset reclustering level2 updates"
```

---

## Task 6: Write final diagnostic figures script

**Files:**
- Create: `pipeline/scripts/11_final_figures.py`
- Input object: `results/12_subset_recluster/adata_annotated_final.h5ad` (output of Task 5, includes updated B cell and erythroid level2)

- [ ] **Step 1: Create the script**

```python
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
```

- [ ] **Step 2: Verify the script compiles**

```bash
pixi run python -m py_compile pipeline/scripts/11_final_figures.py && echo "OK"
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add pipeline/scripts/11_final_figures.py
git commit -m "feat: final figures script for r1.0 annotation (UMAPs, dotplots, attestation)"
```

---

## Task 7: Run final figures and review

**Files:**
- Output: `results/13_final_figures/`

- [ ] **Step 1: Run**

```bash
pixi run python pipeline/scripts/11_final_figures.py \
    --input results/12_subset_recluster/adata_annotated_final.h5ad \
    --out results/13_final_figures
```

Expected runtime: ~3–5 minutes. Expected: `Done. Outputs in results/12_final_figures`

- [ ] **Step 2: Verify key outputs**

```bash
ls results/12_final_figures/ | sort
```

Must include: `umap_manual_level1.png`, `umap_manual_level2.png`, `dotplot_core_markers_by_level1.png`, `dotplot_neutrophil_maturation_by_cluster.png`, `level1_level2_cell_counts.csv`, `ANNOTATION_ATTESTATION.md`

- [ ] **Step 2b: Confirm subset reclustering outcomes are reflected**

```bash
pixi run python -c "
import pandas as pd
df = pd.read_csv('results/13_final_figures/level1_level2_cell_counts.csv')
print(df[df['manual_level1'].isin(['B_cell','erythroid'])].to_string(index=False))
"
```

B_cell should show one of: `{Pro_Pre_B, Mature_B}` (if split was supported) or `{B_cell_broad}`. Erythroid should show `{MEP_Erythroid, Erythroblast}` or `{Erythroid_broad}`.

- [ ] **Step 3: Spot-check the core dotplot**

Open `results/12_final_figures/dotplot_core_markers_by_level1.png`. Verify:
- `granulocyte_neutrophil`: Retnlg/S100a8/S100a9 high, Hlf/Mecom low
- `myeloid_progenitor`: Elane/Mpo/Prtn3 high, Cxcr2/Ccl6 low
- `progenitor`: Hlf/Mecom/Kit high, all myeloid markers low
- `B_cell`: Cd79a/Pax5 high, no myeloid markers
- `DC cluster 20 (pDC)`: Siglech/Bst2 specific

If any compartment shows unexpected cross-reactivity, update the corresponding rows in `pipeline/config/manual_annotation_full_map_leiden_r1.0.tsv` and rerun from Task 3.

- [ ] **Step 4: Review cell count table**

```bash
pixi run python -c "
import pandas as pd
df = pd.read_csv('results/12_final_figures/level1_level2_cell_counts.csv')
print(df.sort_values('n_cells', ascending=False).to_string(index=False))
"
```

- [ ] **Step 5: Commit outputs**

```bash
git add results/13_final_figures/*.csv results/13_final_figures/*.md
git commit -m "results: final annotation figures and attestation (leiden_r1.0 + subset recluster)"
```

---

## Self-Review Checklist

**Spec coverage:**
- [x] Resolution choice justified with concrete cell counts (cluster 35 = 23 cells)
- [x] All 23 leiden_r1.0 clusters annotated at level1
- [x] All 23 clusters annotated at level2 or explicitly flagged as unresolvable
- [x] Sorted gate used as primary biological constraint
- [x] Limitations formally attested (T cell, NK, B cell subtype, HSC/MPP, basophil progenitor)
- [x] Script 10 patched to be resolution-agnostic (reads level2 from map)
- [x] Final figures include neutrophil maturation gradient (the dominant biological signal)

**Placeholder scan:** No TBD, TODO, or "similar to above".

**Type consistency:** `manual_level2` added to `required` set in `read_annotation_map`, mapped via `clusters.map(keyed["manual_level2"])` in `apply_map`, added to `row` dict in `write_cluster_summary`, added to UMAP loop in `plot_umaps`. No naming mismatches.

---

## Summary Table: Level2 Calls by Compartment

| Level1 | Level2 | Confidence | Note |
|--------|--------|------------|------|
| granulocyte_neutrophil | GMP_neutrophil_primed | high | clusters 2, 5 |
| granulocyte_neutrophil | Immature_neutrophil | high | clusters 4, 9, 14, 16 |
| granulocyte_neutrophil | Mature_neutrophil | high/medium | clusters 0, 1, 11, 21 |
| cycling_myeloid | Cycling_myeloid | high | clusters 7, 8 |
| myeloid_progenitor | GMP_neutrophil_primed | high | clusters 2, 5 |
| myeloid_progenitor | Progenitor_ambiguous | low | cluster 6 (heterogeneous CMP/MEP/Bas) |
| progenitor | Progenitor_ambiguous | low | clusters 3, 10 (high mito) |
| monocyte | Classical_monocyte | high | cluster 12 |
| B_cell | B_cell_broad (**subtype unresolvable**) | medium | cluster 15 |
| T_cell | T_cell_broad (**subtype unresolvable**) | high | cluster 18 |
| DC | DC_broad | medium | cluster 17 |
| DC | pDC | high | cluster 20 |
| erythroid | Erythroid_broad (**MEP vs blast unresolvable**) | medium | cluster 13 |
| plasma_cell | Plasma_cell | high | cluster 22 |
| ambiguous_myeloid | Ambiguous (**lineage unresolvable**) | low | cluster 19 |
