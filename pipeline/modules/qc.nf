process QC {
    tag "qc"
    publishDir "${params.outdir}/01_qc", mode: 'copy'
    memory '16 GB'
    cpus 4

    input:
    path h5_file

    output:
    path "adata_qc.h5ad", emit: h5ad
    path "*.png"
    path "cells_per_sample.csv"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/01_qc.py \
        --h5 ${h5_file} \
        --out . \
        --mito_nmads ${params.mito_nmads} \
        --count_nmads ${params.count_nmads}
    """
}
