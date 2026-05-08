process APPLY_MANUAL_ANNOTATION {
    tag "apply_manual_annotation"
    publishDir "${params.outdir}/10_manual_level1", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path h5ad
    path annotation_map

    output:
    path "adata_manual_level1.h5ad", emit: h5ad
    path "*.csv"
    path "*.png"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/10_apply_manual_level1_annotation.py \
        --input ${h5ad} \
        --map ${annotation_map} \
        --out . \
        --cluster-key leiden_r${params.leiden_resolution_manual} \
        --annotation-version ${params.annotation_version}
    """
}
