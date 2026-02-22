process PBMM2_ALIGN {
    tag "$meta.id"
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.14.99--h9ee0642_0':
        'community.wave.seqera.io/library/pbmm2_samtools:c5ded0d11146cfcc' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    val(reset_bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip $fasta -c > ref.fa
    pbmm2 \\
        align \\
        $args \\
        ref.fa \\
        $bam \\
        ${prefix}.bam \\
        --num-threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$(pbmm2 --version |& sed '1!d ; s/pbmm2 //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$(pbmm2 --version |& sed '1!d ; s/pbmm2 //')
    END_VERSIONS
    """
}
