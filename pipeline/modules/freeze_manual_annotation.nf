process FREEZE_MANUAL_ANNOTATION {
    tag "freeze_manual_annotation"
    publishDir "${params.outdir}/10_manual_annotation_frozen", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path h5ad
    path frozen_csv

    output:
    path "adata_progenitor_annotated.h5ad", emit: h5ad
    path "*.csv"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/10b_freeze_manual_annotation.py \
        --input ${h5ad} \
        --frozen ${frozen_csv} \
        --out . \
        --annotation-version ${params.annotation_version}
    """
}
