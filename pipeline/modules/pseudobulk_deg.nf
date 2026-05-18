process PSEUDOBULK_DEG {
    tag "pseudobulk_deg_${level}"
    publishDir "${params.outdir}/19_pseudobulk_deg/${level}", mode: 'copy',
        pattern: 'deg_outputs/*',
        saveAs: { fname -> fname.replaceFirst('^deg_outputs/', '') }
    memory '32 GB'
    cpus 4

    input:
    tuple path(h5ad), val(level)

    output:
    tuple val(level), path("deg_outputs"), emit: deg_dir
    path "deg_outputs/pseudobulk_samples_per_celltype.csv"

    script:
    """
    mkdir -p deg_outputs
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/19_pseudobulk_deg.py \\
        --input ${h5ad} \\
        --out deg_outputs \\
        --level ${level} \\
        --min-cells ${params.pseudobulk_min_cells}
    """
}
