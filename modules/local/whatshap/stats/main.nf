process WHATSHAP_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::whatshap=2.2 bioconda::tabix=1.11"
    container "docker.io/fellen31/whatshap-tabix:2.2"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.stats.tsv") , emit: stats
    tuple val(meta), path("*.blocks.tsv"), emit: blocks
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    whatshap stats \\
        $args \\
        --sample ${meta.id} \\
        --tsv ${prefix}.stats.tsv \\
        --block-list ${prefix}.blocks.tsv \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats.tsv.gz
    touch ${prefix}.blocks.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
    END_VERSIONS
    """
}
