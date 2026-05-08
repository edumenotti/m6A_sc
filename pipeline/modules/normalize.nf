process NORMALIZE {
    tag "normalize"
    publishDir "${params.outdir}/03_normalize", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "adata_normalized.h5ad", emit: h5ad
    path "*.png"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/03_normalize.py \
        --input ${h5ad} \
        --out . \
        --n_hvgs ${params.n_hvgs}
    """
}
