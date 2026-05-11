process CELLCHAT {
    tag "cellchat"
    publishDir "${params.outdir}/16_cellchat", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "cellchat_*.rds"
    path "cellchat_interaction_counts.csv"
    path "cellchat_differential_interactions.csv", optional: true
    path "*.png"

    script:
    """
    pixi run -m ${params.pixi_manifest} Rscript ${projectDir}/scripts/16_cellchat.R \
        --input ${h5ad} \
        --out . \
        --organism ${params.cellchat_organism}
    """
}
