process EXPORT_FOR_R {
    tag "export_for_r"
    publishDir "${params.outdir}/14_progenitor_annotated/r_export", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path h5ad

    output:
    path "r_export", emit: bundle

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/export_h5ad_for_r.py \
        --input ${h5ad} \
        --out r_export \
        --layer counts
    """
}

process CELLCHAT {
    tag "cellchat"
    publishDir "${params.outdir}/16_cellchat", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path bundle

    output:
    path "manual_level1/**", optional: true
    path "manual_level2/**", optional: true

    script:
    """
    pixi run -m ${params.pixi_manifest} Rscript ${projectDir}/scripts/16_cellchat.R \
        --input-dir ${bundle} \
        --out . \
        --organism ${params.cellchat_organism}
    """
}
