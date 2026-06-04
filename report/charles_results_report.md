% m6A perturbation in mouse bone marrow — preliminary single-cell figures
% Charles project — internal, exploratory
% Eduardo Menotti · 2026-06-04

---

**Preliminary analyses — internal discussion only.** Mouse bone-marrow scRNA-seq, ~42,000 cells, 2 × 2 design: genotype (WT, Mutant) × treatment (DMSO vehicle, STM2457 METTL3 inhibitor), two biological donors (paired blocking factor; WT and Mutant separated *in silico*). Cell types use the canonical `cell_type` annotation (HemaScribe lineage calls + manual marker validation, 32 types). Composition was modelled with scCODA, differential expression with pseudobulk PyDESeq2 (`~ donor + genotype + treatment`), and pathway activity with decoupler ULM on the DESeq2 statistic against MSigDB Hallmark/Reactome/KEGG. The genotype contrast is Mutant relative to WT; the treatment contrast is STM relative to DMSO. Figures are descriptive; see caveats at the end.

---

![**Figure 1. Cell-type landscape.** UMAP of ~42,000 bone-marrow cells coloured by the canonical `cell_type` annotation (32 types). Neutrophils are the dominant population (~42%); the stem/progenitor compartment (HSC, STHSC, MPP, GMP) and the lymphoid and erythroid compartments are resolved separately.](figures/fig1_umap_celltypes.png)

![**Figure 2. Cell-type composition.** (A) Per-sample proportion of each cell type on a log scale, coloured by condition (two donors per condition). (B) Mutant-vs-WT log2 fold-change per cell type, pooled over treatment; filled points are flagged credible by scCODA in the genotype contrast, hollow points are not. The stem/progenitor populations HSC and STHSC are lower in Mutant; Monocyte is higher in Mutant. The largest lymphoid fold-changes (B cell, T cell, Pro/Pre-B) coincide with populations captured in only a subset of samples and partly reflect donor/capture differences rather than genotype.](figures/fig2_composition.png)

![**Figure 3. Pathway activity, STM versus DMSO.** Top-40 pathway activities (cell types × pathways) from decoupler ULM. Negative scores indicate lower activity under STM. Translation, ribosome, and rRNA-processing programmes, together with MYC/E2F proliferation programmes, are coordinately lower under STM across progenitor populations.](figures/fig3_pathway_treatment.png)

![**Figure 4. Pathway activity, Mutant versus WT.** Top-40 pathway activities. Negative scores indicate lower activity in Mutant. In HSCs the type-I interferon response is lower in Mutant; the cell-cycle/E2F/G2M direction in HSC is sensitive to how HSC is defined (it inverts between the restrictive HemaScribe and broader manual gates) and is shown for completeness only.](figures/fig4_pathway_genotype.png)

![**Figure 5. Differential expression, MPP4 — Mutant versus WT.** Pseudobulk PyDESeq2 volcano (positive log2FC = higher in Mutant). Genes higher in Mutant include *Irf8*, *Ccl4*, *Il12a*, and *Dntt*; this was the largest genotype DEG response (48 genes at padj < 0.05, |log2FC| > 1).](figures/fig5_volcano_mpp4_genotype.png)

![**Figure 6. Differential expression, cMoP — STM versus DMSO.** Pseudobulk PyDESeq2 volcano (positive log2FC = higher under STM). Genes higher under STM include *Egr1*, *Gbp2*, *Cd74*, and *Mcpt8*.](figures/fig6_volcano_cmop_treatment.png)

---

## Caveats

- Exploratory, hypothesis-generating analyses from a compact pilot dataset — intended to surface directions for follow-up, not to be confirmatory.
- With two donors per group, effect sizes (especially for the treatment contrast) are best read qualitatively rather than as precise estimates.
- The lymphoid compartment is unevenly represented across samples, so lymphoid composition changes are not interpreted here.
- Results for the rarest populations (basophils, megakaryocytes) rest on small numbers and are tentative; scCODA credibility indicates direction, not a calibrated effect size.

*Reproducible from `results/18_composition/final_annotation`, `results/19_pseudobulk_deg/final_annotation`, and `results/20_pathway_activity/final_annotation`. Code: `pipeline/scripts/18_composition_sccoda.py`, `19_pseudobulk_deg.py`, `20_pathway_activity.py`; figures `report/make_fig1_umap.py`, `report/make_fig2_composition.py`. Canonical annotation: `21z_finalize_annotation.py`, `22_clean_annotation_columns.py`. Branch `feat/composition-pseudobulk`.*
