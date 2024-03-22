process TRGT {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/pacbio/trgt:0.4.0"

    input:
    tuple val(meta), path(bam), path(bai), val(sex)
    tuple val(meta2), path(fasta)
    path(repeats)

    output:
    tuple val(meta), path("${meta.id}.spanning.bam"), emit: bam
    tuple val(meta), path("${meta.id}.vcf.gz")      , emit: vcf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    if (sex == '1') {
        karyotype = "XY"
    } else if (sex == '2') {
        karyotype = "XX"
    }

    """
    trgt \\
        $args \\
        --genome ${fasta} \\
        --karyotype ${karyotype} \\
        --repeats ${repeats} \\
        --reads ${bam} \\
        --threads ${task.cpus} \\
        --output-prefix ${meta.id}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trgt: \$(echo \$(trgt -V) | sed 's/TRGT //' )
    END_VERSIONS
    """
}

