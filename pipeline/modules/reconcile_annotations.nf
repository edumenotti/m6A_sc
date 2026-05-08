process RECONCILE_ANNOTATIONS {
    tag "reconcile_annotations"
    publishDir "${params.outdir}/08_reconcile_annotations", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path hspc_h5ad
    path mature_h5ad

    output:
    path "adata_reconciled.h5ad", emit: h5ad
    path "*.csv"
    path "*.json"
    path "*.png"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/08_reconcile_annotations.py \
        --hspc ${hspc_h5ad} \
        --mature ${mature_h5ad} \
        --out . \
        --cluster-key leiden_r${params.leiden_resolution} \
        --hspc-label-regex '${params.reconcile_hspc_label_regex}' \
        --min-confidence ${params.reconcile_min_confidence} \
        --cluster-majority-min-fraction ${params.reconcile_cluster_majority_min_fraction}
    """
}
