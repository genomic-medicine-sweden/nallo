process ECHTVAR_ANNO {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/fellen31/echtvar:0.2.2"

    input:
    tuple val(meta),  path(vcf)
    path(databases)

    output:
    tuple val(meta), path("*.bcf.gz"), emit: bcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def input = databases.collectMany { file -> ["-e", file] }.join(" ")
    """
    echtvar \\
        anno \\
        ${args} \\
        ${input} \\
        ${vcf} \\
        ${prefix}.bcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echtvar: \$(echo \$(echtvar -V) | sed 's/echtvar //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echtvar: \$(echo \$(echtvar -V) | sed 's/echtvar //' )
    END_VERSIONS
    """
}
