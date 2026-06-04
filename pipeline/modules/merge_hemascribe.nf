process MERGE_HEMASCRIBE {
    tag "merge_hemascribe"
    publishDir "${params.outdir}/21_hemascribe", mode: 'copy',
        pattern: 'audit/*',
        saveAs: { fname -> fname.replaceFirst('^audit/', '') }
    memory '32 GB'
    cpus 2

    input:
    tuple path(h5ad), path(labels)

    output:
    path "adata_hemascribe_full.h5ad", emit: h5ad
    path "audit/*", optional: true

    script:
    """
    mkdir -p audit
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/21c_merge_hemascribe_labels.py \
        --adata ${h5ad} \
        --labels ${labels} \
        --out adata_hemascribe_full.h5ad \
        --audit-dir audit
    """

    stub:
    """
    mkdir -p audit
    touch adata_hemascribe_full.h5ad audit/confusion_stub.csv
    """
}
