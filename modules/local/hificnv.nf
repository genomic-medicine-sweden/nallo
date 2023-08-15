process HIFICNV {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::hificnv=0.1.6b"
    container "quay.io/biocontainers/hificnv:0.1.6b--h9ee0642_0"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    // Take care of this later
    tuple val(meta3), path(maf_vcf_gz)
    path(exclude_bed)
    path(expected_cn)

    output:
    tuple val(meta), path("*.vcf.gz")  , emit: vcf
    tuple val(meta), path("*.depth.bw"), emit: depth_bw
    tuple val(meta), path("*.maf.bw")  , emit: maf_bw, optional: true
    tuple val(meta), path("*.bedgraph"), emit: cnval
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    // TODO: add exclude, expected
    """
    hificnv \\
        $args \\
        --bam ${bam} \\
        --ref ${fasta} \\
        --threads $task.cpus \\
        --output-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(echo \$(cramino -V) | sed 's/cramino //' )
    END_VERSIONS
    """
}

