process CLUSTER {
    tag "cluster"
    publishDir "${params.outdir}/05_cluster", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "adata_clustered.h5ad", emit: h5ad
    path "*.png"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/05_cluster.py \
        --input ${h5ad} \
        --out . \
        --resolution ${params.leiden_resolution}
    """
}
