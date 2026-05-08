process MARKERS {
    tag "markers"
    publishDir "${params.outdir}/07_markers", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "*.png"
    path "*.csv"

    script:
    """
    pixi run -m ${projectDir}/../pixi.toml python ${projectDir}/scripts/07_markers.py \
        --input ${h5ad} \
        --out . \
        --cluster_key leiden_r${params.leiden_resolution} \
        --marker-db ${projectDir}/config/cell_type_markers.tsv \
        --watchlist-clusters 10,12
    """
}
