process MACROPHAGE_STATES {
    tag "macrophage_states"
    publishDir "${params.outdir}/15_macrophage_states", mode: 'copy'
    memory '16 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "monocyte_state_scores.csv"
    path "monocyte_score_summary.csv"
    path "*.png"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/15_macrophage_states.py \
        --input ${h5ad} \
        --out .
    """
}
