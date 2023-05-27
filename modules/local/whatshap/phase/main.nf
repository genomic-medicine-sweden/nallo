process WHATSHAP_PHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::whatshap=1.7 bioconda::tabix=1.11"
    container "docker.io/fellen31/whatshap-tabix:latest"

    input:
    tuple val(meta),  path(vcf)
    tuple val(meta2), path(bam), path(bai)
    tuple val(meta3), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path(".command.err"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    whatshap phase \\
        $args \\
        -o ${vcf}.phased.vcf \\
        --reference $fasta \\
        ${vcf} \\
        ${bam}

    bgzip \\
        -@ $task.cpus \\
        ${vcf}.phased.vcf

    tabix \\
        -p vcf \\
        ${vcf}.phased.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
        bgzip: \$( bgzip --version | head -n 1 | sed 's/bgzip (htslib) //g')
        tabix: \$( tabix --version | head -n 1 | sed 's/tabix (htslib) //g')
    END_VERSIONS
    """
}
