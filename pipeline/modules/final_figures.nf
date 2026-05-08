process FINAL_FIGURES {
    tag "final_figures"
    publishDir "${params.outdir}/11_final_figures", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path h5ad

    output:
    path "*.png"
    path "*.csv"
    path "*.md"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/11_final_figures.py \
        --input ${h5ad} \
        --out .
    """
}
