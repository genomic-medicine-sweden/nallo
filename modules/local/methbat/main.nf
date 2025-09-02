process METHBAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methbat:0.1.0--pyhdfd78af_0':
        'biocontainers/methbat:0.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)
    path(bed_regions), stageAs: "regions.bed"

    output:
    tuple val(meta), path("*.tsv"),         emit: methylation_calls
    tuple val(meta), path("*.bed"),         emit: bed,              optional: true
    tuple val(meta), path("*.summary.txt"), emit: summary
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_arg = bed_regions ? "--regions ${bed_regions}" : ""
    
    """
    methbat \\
        --bam ${bam} \\
        --ref ${reference} \\
        --output-prefix ${prefix} \\
        ${regions_arg} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methbat: \$(methbat --version 2>&1 | grep -o 'MethBat [0-9.]*' | cut -d' ' -f2 || echo "0.1.0")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.methylation_calls.tsv
    touch ${prefix}.regions.bed
    touch ${prefix}.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methbat: 0.1.0
    END_VERSIONS
    """
}