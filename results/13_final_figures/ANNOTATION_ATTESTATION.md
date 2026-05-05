# Annotation Attestation — Charles Mouse BM scRNA-seq
# leiden_r1.0, 42390 cells, 11 level1 categories

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
