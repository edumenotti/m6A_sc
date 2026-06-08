# m6A_sc — bone marrow scRNA-seq (competitive transplant)

Single-cell RNA-seq of a paired competitive bone-marrow transplant
(WT vs Mutant, ± STM): from raw counts to an annotated object plus downstream
analyses (scCODA composition, PyDESeq2 pseudobulk DEG, HemaScribe validation).

The whole analysis is a [Nextflow](https://www.nextflow.io/) pipeline
(`pipeline/main.nf`) running inside one Apptainer/Singularity container
(`charles-scrna.def`) with a pinned [pixi](https://pixi.sh) env (`pixi.toml`).

## Run

```bash
# build the container once
apptainer build charles-scrna.sif charles-scrna.def

# run the full pipeline (Yale Bouchet HPC: Slurm + Apptainer)
nextflow run pipeline/main.nf -profile bouchet -params-file pipeline/params.yaml
```

Other profiles: `local` (host pixi env), `slurm`, `apptainer`.
Parameters: `pipeline/params.yaml`.

## Output

- Annotated object: `results/14_progenitor_annotated/adata_hemascribe.h5ad`
- Results write-up: `report/charles_results_report.html`
