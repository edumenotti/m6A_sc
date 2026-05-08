process MANUAL_MARKER_REVIEW {
    tag "manual_marker_review"
    publishDir "${params.outdir}/09_manual_marker_review", mode: 'copy'
    memory '24 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "*.png"
    path "*.csv"
    path "*.tsv"
    path "*.md"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/09_manual_marker_review.py \
        --input ${h5ad} \
        --out . \
        --marker-db ${projectDir}/config/manual_annotation_markers_skull_immune.tsv \
        --cluster-key leiden_r${params.leiden_resolution_manual}
    """
}
