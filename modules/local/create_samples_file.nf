process CREATE_SAMPLES_FILE {
    tag "${meta.id}"
    label 'process_single'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("*.txt"), emit: samples
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo SAMPLE "${meta.id}" > ${meta.id}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_samples_file: v1.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_samples_file: v1.0
    END_VERSIONS
    """
}
