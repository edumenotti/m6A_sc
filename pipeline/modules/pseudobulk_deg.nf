process PSEUDOBULK_DEG {
    tag "pseudobulk_deg_${level}"
    publishDir "${params.outdir}/19_pseudobulk_deg/${level}", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    tuple path(h5ad), val(level)

    output:
    path "deg_*"
    path "pseudobulk_samples_per_celltype.csv"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/19_pseudobulk_deg.py \\
        --input ${h5ad} \\
        --out . \\
        --level ${level} \\
        --min-cells ${params.pseudobulk_min_cells}
    """
}
