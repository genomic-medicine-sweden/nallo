process CRAMINO {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::cramino=0.14.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cramino:0.14.5--h5076881_0' :
        'biocontainers/cramino:0.14.5--h5076881_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.txt"),   emit: stats
    tuple val(meta), path("*.arrow"), emit: arrow
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cramino \\
        $args \\
        --threads $task.cpus \\
        --arrow ${bam.baseName}.arrow \\
        ${bam} > ${bam.baseName}.cramino.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: \$(echo \$(cramino -V) | sed 's/cramino //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${bam.baseName}.arrow
    touch ${bam.baseName}.cramino.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: \$(echo \$(cramino -V) | sed 's/cramino //' )
    END_VERSIONS
    """
}

