# Manual marker-aware annotation plan

## Goal

Replace `popV` as the final authority with an auditable manual annotation workflow for the Charles mouse bone marrow scRNA-seq dataset. Keep `popV`/scANVI labels as auxiliary evidence, but assign final labels from marker programs, DE genes, QC overlays, and cluster/subcluster structure.

## Reference workflow

Primary reference: `external/skull_immune/transcriptomics/02_celltype_annotation/`.

The reusable pattern from `skull_immune` is:

1. Annotate coarse clusters first.
2. For each coarse population, create a focused subset.
3. Recompute neighborhood/UMAP/leiden on the subset when a compartment is heterogeneous.
4. Review known marker genes, top DE genes, QC covariates, and sample/condition composition.
5. Assign explicit `level1` and `level2` labels through a manual mapping table.
6. Preserve excluded/low-quality/hybrid decisions as flags or separate review categories, not as ordinary cell-type labels.

## Charles adaptation

Use `results/08_reconcile_annotations/adata_reconciled.h5ad` as the working object because it contains preprocessing, scVI integration, clustering, popV/scANVI evidence, and reconciled labels.

Use `leiden_r1.0` as the initial manual review granularity. It already splits the problematic `leiden_r0.5` clusters well enough for first-pass marker correction.

Keep these columns:

- `final_cell_type`: current reconciled automated label.
- popV/scANVI prediction columns: auxiliary evidence only.
- QC metrics and marker scores: audit evidence.

Add new columns rather than overwriting:

- `manual_level1`
- `manual_level2`
- `manual_annotation_basis`
- `manual_annotation_confidence`
- `manual_review_flag`

## Marker catalog

Candidate marker table: `pipeline/config/manual_annotation_markers_skull_immune.tsv`.

This table merges:

- `skull_immune` dotplot markers from `02e_final_annotations_dotplot.ipynb`.
- Existing Charles markers from `pipeline/config/cell_type_markers.tsv`.
- Charles `leiden_r1.0` DE evidence from `results/08_reconcile_annotations/reconciled_cluster_marker_review_leiden_r1.0.csv`.

The table is intentionally broader than the first-pass labels. Markers must be validated by expression in this dataset before they are used for hard assignment.

## First-pass label decisions to test

Use the `leiden_r1.0` review table to seed an explicit manual map:

- `leiden_r1.0` 0, 1, 4, 9, 11, 14, 16, 21: neutrophil/granulocyte-like, with level2 maturity refined by markers.
- `leiden_r1.0` 12: monocyte-like.
- `leiden_r1.0` 15: B-cell-like.
- `leiden_r1.0` 18: T-cell-like.
- `leiden_r1.0` 22: plasma-cell-like.
- `leiden_r1.0` 3: HSC/MPP-like, but this cluster needs careful QC/ribosomal review before finalizing.

Do not treat `Low_Quality_Ambient`, `Stress_Response`, or `Cycling` as cell-type labels. Use them as flags.

## Implementation phases

1. Score the expanded marker catalog on the reconciled object.
2. Generate a cluster review table for all `leiden_r1.0` clusters, not only the prior watchlist.
3. Draft `pipeline/config/manual_annotation_map_leiden_r1.0.tsv`.
4. Write a script that applies the manual map and creates `manual_level1`/`manual_level2`/audit columns.
5. Generate final UMAPs, dotplots, cluster composition tables, and disagreement tables against popV/scANVI.
6. Review ambiguous compartments with targeted subset reclustering only if the full `r1.0` evidence is insufficient.

## Open questions

- Whether neutrophil-like clusters should be split into pro/pre/immature/mature neutrophils in the first manual pass or kept as broader neutrophil-like labels until subset reclustering.
- Whether lymphoid/Input-enriched clusters represent true biological carryover, expected marrow lymphocytes, or sample-specific contamination.
- Whether erythroid and megakaryocyte compartments need targeted reclustering before final labels.
- Whether a small number of `DC`/pDC clusters remain valid after removing granulocyte/monocyte overcalls.
