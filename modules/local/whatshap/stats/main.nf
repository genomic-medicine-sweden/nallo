process WHATSHAP_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::whatshap=1.7 bioconda::tabix=1.11"
    container "docker.io/fellen31/whatshap-tabix:latest"

    input:
    tuple val(meta),  path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.stats.tsv.gz"), emit: stats
    tuple val(meta), path("*.blocks.tsv")  , emit: blocks
    tuple val(meta), path(".command.err")  , emit: err
    tuple val(meta), path(".command.log")  , emit: log
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    whatshap stats \\
        $args \\
        --sample ${meta.id} \\
        --tsv ${vcf}.stats.tsv.gz \\
        --block-list ${vcf}.blocks.tsv \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
    END_VERSIONS
    """
}
