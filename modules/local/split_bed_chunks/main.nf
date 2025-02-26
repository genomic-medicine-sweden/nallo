process SPLIT_BED_CHUNKS {
    tag "$meta.id"

    container "quay.io/biocontainers/pandas:1.5.2"

    input:
    tuple val(meta), path(bed)
    val(chunk_size)

    output:
    path("*.bed")       , emit: split_beds
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // Chunk size needs to be > 0. Where to specify a needed argument? params?
    // Will not output more regions than in file, even if chunk_size > number of regions
    script:
    def VERSION = "1.0"
    """
    split_bed_chunks.py ${bed} $chunk_size

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split_bed_chunks: \$(echo "$VERSION" )
    END_VERSIONS
    """
}
