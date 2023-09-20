process ECHTVAR_ANNO {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/fellen31/echtvar:latest"

    input:
    tuple val(meta),  path(vcf)
    path(databases)

    output:
    tuple val(meta), path("${meta.id}.bcf.gz"), emit: bcf
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    modifiedList = []
    for (element in databases) {
        modifiedList.add("-e")
        modifiedList.add(element)
    }
    println databases
    """
    echtvar anno ${modifiedList.join(" ")} ${vcf} ${meta.id}.bcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echtvar: \$(echo \$(echtvar -V) | sed 's/echtvar //' )
    END_VERSIONS
    """
}

