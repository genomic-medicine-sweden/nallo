process BEDMETHYL_CONCAT {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5':
        'biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path( "${output_file}" )   , emit: sorted
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    suffix          = task.ext.suffix   ?: '.bed'
    def input_list  = input.sort{it.toString()}.join(" ")
    output_file     = "${prefix}_${meta.group}.${suffix}"
    if ("$input" == "$output_file") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    sort \\
        --parallel ${task.cpus} \\
        ${args} \\
        ${input_list} > ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version |& sed '1!d ; s/sort (GNU coreutils) //')
    END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    suffix          = task.ext.suffix   ?: "${input.extension}"
    output_file     = "${prefix}.${suffix}"
    if ("$input" == "$output_file") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version |& sed '1!d ; s/sort (GNU coreutils) //')
    END_VERSIONS
    """
}
