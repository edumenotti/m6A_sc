process PATHWAY_ACTIVITY {
    tag "pathway_activity_${level}"
    publishDir "${params.outdir}/20_pathway_activity/${level}", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    tuple path(deg_dir), val(level)

    output:
    path "pathway_*"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/20_pathway_activity.py \\
        --deg-dir ${deg_dir} \\
        --out . \\
        --organism ${params.pathway_organism}
    """
}
