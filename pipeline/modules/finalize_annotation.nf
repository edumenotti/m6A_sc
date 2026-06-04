process FINALIZE_ANNOTATION {
    tag "finalize_annotation"
    publishDir "${params.outdir}/21_hemascribe", mode: 'copy',
        pattern: 'audit/*',
        saveAs: { fname -> fname.replaceFirst('^audit/', '') }
    memory '32 GB'
    cpus 2

    input:
    path h5ad

    output:
    path "adata_hemascribe_finalized.h5ad", emit: h5ad
    path "audit/*", optional: true

    script:
    """
    mkdir -p audit
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/21z_finalize_annotation.py \
        --input ${h5ad} \
        --out adata_hemascribe_finalized.h5ad \
        --audit-dir audit
    """

    stub:
    """
    mkdir -p audit
    touch adata_hemascribe_finalized.h5ad audit/final_annotation_summary.csv
    """
}
