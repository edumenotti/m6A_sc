process EXPORT_FOR_HEMASCRIBE {
    tag "export_for_hemascribe"
    publishDir "${params.outdir}/21_hemascribe/export", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    path h5ad

    output:
    path "hemascribe_export", emit: export_dir

    script:
    """
    mkdir -p hemascribe_export
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/21a_export_for_hemascribe.py \
        --input ${h5ad} \
        --out hemascribe_export \
        --layer counts
    """

    stub:
    """
    mkdir -p hemascribe_export
    touch hemascribe_export/counts.mtx hemascribe_export/barcodes.tsv.gz \
          hemascribe_export/features.tsv.gz hemascribe_export/obs.csv
    """
}
