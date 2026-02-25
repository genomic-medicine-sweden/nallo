

process PARAPHRASE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'docker.io/fellen31/paraphrase:0.1.0' }"

    input:
    tuple val(meta), path(jsons), val(samples)
    tuple val(meta2), path(yaml)
    val(output_format)

    output:
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val(meta), path("*.tsv"), emit: tsv, optional: true

    tuple val("${task.process}"), val('paraphrase'), eval("paraphrase --version"), topic: versions, emit: versions_paraphrase

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rules = yaml ? "--rules $yaml" : ''
    """
    paraphrase \
        $args \
        --input ${jsons.join(' --input ')} \
        --sample ${samples.join(' --sample ')} \
        $rules \
        --output-format $output_format \
        > ${prefix}.${output_format}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    """
}
