process CLEAN_COLUMNS {
    tag "clean_columns"
    publishDir "${params.outdir}/14_progenitor_annotated", mode: 'copy'
    memory '32 GB'
    cpus 2

    input:
    path h5ad

    output:
    path "adata_hemascribe.h5ad", emit: h5ad

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/22_clean_annotation_columns.py \
        --input ${h5ad} \
        --out adata_hemascribe.h5ad
    """

    stub:
    """
    touch adata_hemascribe.h5ad
    """
}
