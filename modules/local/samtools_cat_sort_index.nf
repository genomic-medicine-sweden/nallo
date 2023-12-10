process SAMTOOLS_CAT_SORT_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta),  path(input_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${prefix}.bam") , optional:true, emit: bam
    tuple val(meta), path("${prefix}.bam.bai") , optional:true, emit: bai
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai") , optional:true, emit: bam_bai
    path  "versions.yml"                                  , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def args2 = task.ext.args2   ?: ''
    def args3 = task.ext.args3   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    """
    samtools \\
        cat \\
        --threads ${task.cpus-1} \\
        $args \\
        $input_files |\\
    samtools \\
        sort \\
        --threads ${task.cpus-1} \\
        $args2 \\
        -O BAM \\
        -o ${prefix}.${file_type} \\

    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        ${prefix}.${file_type} \\
        $args3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    """
    touch ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
