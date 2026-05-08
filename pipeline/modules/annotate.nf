process ANNOTATE {
    tag "annotate"
    publishDir "${params.outdir}/06_annotate", mode: 'copy'
    memory '32 GB'
    cpus 4
    accelerator 1

    input:
    path h5ad
    path ref_h5ad

    output:
    path "adata_annotated.h5ad", emit: h5ad
    path "*.png"
    path "*.csv"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml -e popv python ${projectDir}/scripts/06_annotate.py \
        --input ${h5ad} \
        --ref ${ref_h5ad} \
        --out . \
        --ref-labels-key ${params.ref_labels_key} \
        --ref-batch-key ${params.ref_batch_key} \
        --query-batch-key sample_id \
        --cluster-key leiden_r${params.leiden_resolution} \
        --popv-methods ${params.popv_methods} \
        --popv-mode ${params.popv_mode} \
        --popv-hvg ${params.popv_hvg} \
        --popv-samples-per-label ${params.popv_samples_per_label} \
        --scanvi-n-latent ${params.n_latent}
    """
}
