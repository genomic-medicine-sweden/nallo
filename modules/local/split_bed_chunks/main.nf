process SPLIT_BED_CHUNKS {
    tag "$intervals"

    container "quay.io/biocontainers/pandas:1.5.2"

    input:
    tuple val(meta), path(bed)
    val(chunk_size)

    output:
    path("*.bed")       , emit: split_beds
    //path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // Chunk size needs to be > 0. Where to specify a needed argument? params?
    script:
    """
    split_bed_chunks.py ${bed} $chunk_size
    """
}
