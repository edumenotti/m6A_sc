process COMPOSITION_SCCODA {
    tag "composition_${level}"
    publishDir "${params.outdir}/18_composition/${level}", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    tuple path(h5ad), val(level)

    output:
    path "sccoda_*"

    script:
    """
    export CUDA_VISIBLE_DEVICES=""
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/18_composition_sccoda.py \\
        --input ${h5ad} \\
        --out . \\
        --level ${level}
    """
}
