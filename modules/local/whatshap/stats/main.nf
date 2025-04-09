process WHATSHAP_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::whatshap=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:2.3--py38h2494328_0' :
        'quay.io/biocontainers/whatshap:2.3--py38h2494328_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.stats.tsv")        , emit: stats
    tuple val(meta), path("*.blocks.tsv.gz")    , emit: blocks
    tuple val(meta), path("*.blocks.tsv.gz.tbi"), emit: blocks_index
    path "versions.yml"                         , emit: versions

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

    bgzip \\
        -@ $task.cpus \\
        ${prefix}.blocks.tsv

    tabix \\
        ${prefix}.blocks.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
        bgzip: \$( bgzip --version | head -n 1 | sed 's/bgzip (htslib) //g')
        tabix: \$( tabix --version | head -n 1 | sed 's/tabix (htslib) //g')
    END_VERSIONS
    """
    stub:

    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.stats.tsv
    echo | gzip > ${prefix}.blocks.tsv.gz
    touch ${prefix}.blocks.tsv.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
        bgzip: \$( bgzip --version | head -n 1 | sed 's/bgzip (htslib) //g')
        tabix: \$( tabix --version | head -n 1 | sed 's/tabix (htslib) //g')
    END_VERSIONS
    """
}
