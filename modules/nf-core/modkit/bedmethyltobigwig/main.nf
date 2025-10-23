process MODKIT_BEDMETHYLTOBIGWIG {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/fellen31/modkit:v0.5.1-rc1"

    input:
    tuple val(meta),  path(bedmethyl)
    tuple val(meta2), path(chromsizes)
    val modcodes

    output:
    tuple val(meta), path("*.bw"), emit: bw
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mods = modcodes instanceof List ? modcodes.join(',') : modcodes
    """
    modkit bedmethyl tobigwig \\
        $args \\
        --nthreads $task.cpus \\
        --sizes $chromsizes \\
        --mod-codes $mods \\
        $bedmethyl \\
        ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | sed 's/modkit //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | sed 's/modkit //')
    END_VERSIONS
    """
}
