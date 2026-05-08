process DOUBLETS {
    tag "doublets"
    publishDir "${params.outdir}/02_doublets", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "adata_no_doublets.h5ad", emit: h5ad
    path "*.png"
    path "*.csv"

    script:
    def installFlag = params.doubletfinder_install_missing ? "--install-missing" : ""
    def autoPkFlag = params.doubletfinder_auto_pk ? "--auto-pk" : ""
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/02_doublets.py \
        --input ${h5ad} \
        --out . \
        --expected-rate ${params.expected_doublet_rate} \
        --pcs ${params.doubletfinder_pcs} \
        --pk ${params.doubletfinder_pk} \
        --num-cores ${task.cpus} \
        ${installFlag} \
        ${autoPkFlag}
    """
}
