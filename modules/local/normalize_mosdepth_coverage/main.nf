process NORMALIZE_MOSDEPTH_COVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'python:3.11' }"

    input:
    tuple val(meta), path(mosdepth_tsv)

    output:
    tuple val(meta), path("*.tsv"), emit: normalized
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    normalize_mosdepth_coverage.py \
        --input ${mosdepth_tsv} \
        --output ${prefix}.normalized.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        normalize_mosdepth_coverage: \$(normalize_mosdepth_coverage.py --version | sed 's/normalize_mosdepth_coverage //')
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.normalized.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        normalize_mosdepth_coverage: \$(normalize_mosdepth_coverage.py --version | sed 's/normalize_mosdepth_coverage //')
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
