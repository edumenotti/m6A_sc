process APPLY_PROGENITOR_ANNOTATION {
    tag "apply_progenitor_annotation"
    publishDir "${params.outdir}/14_progenitor_annotated", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path h5ad
    path assignments
    path map_tsv

    output:
    path "adata_progenitor_annotated.h5ad", emit: h5ad
    path "final_annotation_counts.csv"
    path "umap_final_annotation.png"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/14_apply_progenitor_annotation.py \\
        --input ${h5ad} \\
        --assignments ${assignments} \\
        --map ${map_tsv} \\
        --out .
    """
}
