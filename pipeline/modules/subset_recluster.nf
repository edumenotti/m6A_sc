process SUBSET_RECLUSTER {
    tag "subset_recluster"
    publishDir "${params.outdir}/12_subset_recluster", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "adata_annotated_final.h5ad", emit: h5ad
    path "*.png"
    path "*.csv"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/12_subset_recluster.py \
        --input ${h5ad} \
        --out .
    """
}
