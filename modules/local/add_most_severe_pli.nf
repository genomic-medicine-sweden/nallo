process ADD_MOST_SEVERE_PLI {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$vcf" == "${prefix}.vcf" ) error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    add_most_severe_pli.py \\
        --file_in ${vcf} \\
        --file_out ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_most_severe_pli: 1.1
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pli.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_most_severe_pli: 1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
