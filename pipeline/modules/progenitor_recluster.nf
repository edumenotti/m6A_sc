process PROGENITOR_RECLUSTER {
    tag "progenitor_recluster"
    publishDir "${params.outdir}/13_progenitor_recluster", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "annotation_template.tsv", emit: template
    path "leiden_subset_assignments.csv", emit: assignments
    path "*.csv"
    path "*.png"
    path "*.tsv"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/13_progenitor_recluster.py \\
        --input ${h5ad} \\
        --out .
    """
}
