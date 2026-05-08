process INTEGRATE {
    tag "integrate"
    publishDir "${params.outdir}/04_integrate", mode: 'copy'
    memory '64 GB'
    cpus 8
    accelerator 1

    input:
    path h5ad

    output:
    path "adata_integrated.h5ad", emit: h5ad
    path "scvi_model/", emit: model
    path "*.png"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/04_integrate.py \
        --input ${h5ad} \
        --out . \
        --n_latent ${params.n_latent} \
        --max_epochs ${params.max_epochs}
    """
}
