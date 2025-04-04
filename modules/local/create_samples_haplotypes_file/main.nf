process CREATE_SAMPLES_HAPLOTYPES_FILE {
    tag "${meta.id}"
    label 'process_single'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path("*.txt"), emit: samples
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$txt" == "${prefix}.txt") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    awk '{print \$1,"${meta.id}_"\$1}' ${txt} > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
