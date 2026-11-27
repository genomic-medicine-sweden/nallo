process SNIFFLES1 {
    tag "${meta.id}"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sniffles:1.0.12--h8b12597_0'
        : 'biocontainers/sniffles:1.0.12--h8b12597_0'}"

    input:
    tuple val(meta), path(input), path(bai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('sniffles'), val('1.0.12'), emit: versions_sniffles1, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sniffles \\
        --mapped_reads ${input} \\
        --threads ${task.cpus} \\
        -v ${prefix}.vcf \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
