# m6A_sc — bone marrow scRNA-seq (competitive transplant)

Single-cell RNA-seq analysis of a **paired competitive bone-marrow transplant**:
WT vs Mutant cells, ± STM treatment. The pipeline goes from raw counts to a
fully annotated AnnData object and the downstream analyses reported to the lab.

## What's in here

- **End-to-end [Nextflow](https://www.nextflow.io/) pipeline** (`pipeline/main.nf`)
  running inside a single **Apptainer/Singularity** container
  (`charles-scrna.def`) with a pinned [pixi](https://pixi.sh) environment
  (`pixi.toml`). One command, reproducible on any HPC — no manual setup.
- **Pipeline stages:** QC → doublet removal → normalization → scVI integration →
  Leiden clustering → (optional reference label transfer) → manual annotation
  (transferred per-barcode for reproducibility) → HemaScribe annotation
  validation → canonical `cell_type` (32 types).
- **Downstream analyses** on the canonical annotation: compositional analysis
  (scCODA), per-celltype pseudobulk DEG (PyDESeq2, donor as replicate), pathway
  activity, and per-cell macrophage-state scoring.
- **Report:** `report/charles_results_report.html` (and `.pdf`) — the results
  write-up with figures.

## How to run

Build the container once:

```bash
apptainer build charles-scrna.sif charles-scrna.def
```

Run the full pipeline (Yale Bouchet HPC — Slurm + Apptainer in one profile):

```bash
nextflow run pipeline/main.nf -profile bouchet -params-file pipeline/params.yaml
```

Other profiles: `standard`/`local` (host pixi env, no container), `slurm`,
`apptainer`. Pipeline parameters live in `pipeline/params.yaml`.

## Final output

The canonical annotated object is written to:

```
results/14_progenitor_annotated/adata_hemascribe.h5ad
```

42,384 cells × 32,285 genes, 32 `cell_type` labels, WT/Mutant genotype per cell.

## Reproducibility note

The full pipeline was re-run on the HPC and reproduces the local result
cell-for-cell on cell-type composition and genotype assignment (differences are
< 0.02% of cells, from one stochastic upstream QC step; all annotation labels are
transferred by barcode and reproduce exactly). Embeddings (scVI/UMAP) and Leiden
cluster IDs are not expected to be bit-identical across hardware — the annotation
is intentionally barcode-anchored so it survives that drift.

## Data notes

Design is n=2 donors with a sort confound and an in-silico genotype split;
cross-condition claims are interpreted with that ceiling in mind. Exploratory
analyses that didn't survive those caveats are archived under
`results/_archive_exploratory/`.
