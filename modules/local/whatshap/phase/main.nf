process WHATSHAP_PHASE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::whatshap=2.2 bioconda::tabix=1.11"
    container "docker.io/fellen31/whatshap-tabix:2.2"

    input:
    tuple val(meta), path(vcf), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_tbi
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    whatshap phase \\
        $args \\
        -o ${prefix}.vcf \\
        --reference $fasta \\
        ${vcf} \\
        ${bam}

    bgzip \\
        -@ $task.cpus \\
        ${prefix}.vcf

    tabix \\
        -p vcf \\
        ${prefix}.vcf.gz

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
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
        bgzip: \$( bgzip --version | head -n 1 | sed 's/bgzip (htslib) //g')
        tabix: \$( tabix --version | head -n 1 | sed 's/tabix (htslib) //g')
    END_VERSIONS
    """
}
