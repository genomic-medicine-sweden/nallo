process ECHTVAR_ENCODE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/fellen31/echtvar:latest"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("${meta.id}.zip"), emit: db
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    """
    cat <<-JSON > ${prefix}.json
    [
        {
            "field": "AF",
            "alias": "COHORT_AF",
            "multiplier": 1000000
        },
        {
            "field": "AC",
            "alias": "COHORT_AC",
            "multiplier": 1000000
        },

    ]
    JSON

    echtvar encode ${prefix}.zip ${prefix}.json ${bcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echtvar: \$(echo \$(echtvar -V) | sed 's/echtvar //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        echtvar: \$(echo \$(echtvar -V) | sed 's/echtvar //' )
    END_VERSIONS
    """

}

