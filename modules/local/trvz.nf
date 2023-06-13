process TRVZ {
    tag "$meta.id"
    label 'process_single'
    
    container "quay.io/pacbio/trgt:0.4.0"

    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(csi), path(single_repeat)
    tuple val(meta2), path(fasta)
    path(repeats)

    output:
    tuple val(meta), path("*.svg"), emit: svg
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    
    """
    trvz \\
        $args \\
        --genome ${fasta} \\
        --repeats ${repeats} \\
        --spanning-reads ${bam} \\
        --vcf ${vcf} \\
        --repeat-id ${single_repeat.baseName} \\
        --image ${meta.id}.${single_repeat.baseName}.svg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trgt: \$(echo \$(trgt -V) | sed 's/TRGT //' )
    END_VERSIONS
    """
}

