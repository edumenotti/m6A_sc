process NICHENET {
    tag "nichenet"
    publishDir "${params.outdir}/17_nichenet", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "nichenet_ligand_activity_*.csv"
    path "*.png"

    script:
    """
    pixi run -m ${params.pixi_manifest} Rscript ${projectDir}/scripts/17_nichenet.R \
        --input ${h5ad} \
        --out . \
        --organism ${params.nichenet_organism}
    """
}
