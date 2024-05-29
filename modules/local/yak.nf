process YAK {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::yak=0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yak:0.1--h7132678_2':
        'biocontainers/yak:0.1--h7132678_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.yak"), emit: yak
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    yak \\
        count \\
        -k31 \\
        -b37 \\
        $args \\
        -t $task.cpus \\
        -o ${prefix}.yak \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yak: \$(echo yak version))
    END_VERSIONS
    """
}
