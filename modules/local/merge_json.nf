process MERGE_JSON {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(json_files)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_json_sample_files_into_family.py \\
        --files_in ${json_files} \\
        --file_out ${prefix}_merged.json \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_json_sample_files_into_case: 1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_json_sample_files_into_family: 1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
