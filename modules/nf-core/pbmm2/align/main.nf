process PBMM2_ALIGN {
    tag "$meta.id"
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b753fe18151c603551fb74c541421bf1ceb4249ea86de9f2ede6120ccf6d8f7/data':
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
    samtools reset -@ ${task.cpus} $bam -O BAM > reset.bam
    pbmm2 \\
        align \\
        $args \\
        ref.fa \\
        reset.bam \\
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
