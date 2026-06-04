process HEMASCRIBE {
    tag "hemascribe"
    publishDir "${params.outdir}/21_hemascribe", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path export_dir

    output:
    path "hemascribe_labels.csv", emit: labels

    script:
    """
    pixi run -m ${params.pixi_manifest} Rscript ${projectDir}/scripts/21b_hemascribe.R \
        --indir ${export_dir} \
        --out hemascribe_labels.csv
    """

    stub:
    """
    echo "cell_id,broad.annot,fine.annot,hematopoietic.score" > hemascribe_labels.csv
    """
}
