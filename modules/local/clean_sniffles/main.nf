process CLEAN_SNIFFLES {
    tag "${meta.id}"

    conda "conda-forge::python=3.8.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.8.3'
        : 'biocontainers/python:3.8.3'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0"
    """
    clean_sniffles.py ${vcf} > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clean_sniffles: \$(echo "${VERSION}" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clean_sniffles: \$(echo "${VERSION}" )
    END_VERSIONS
    """
}
