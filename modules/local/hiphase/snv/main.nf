process HIPHASE_SNV {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/biocontainers/hiphase:0.10.0--h9ee0642_0"

    input:
    // bgzipped, indexed inputs, sample name needs to be match + Sniffles fill REF tags
    tuple val(meta), path(snp_vcf), path(snp_csi),path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hiphase \\
        $args \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        --bam ${bam} \\
        --vcf ${snp_vcf} \\
        --output-vcf ${snp_vcf.baseName}.phased.vcf.gz \\
        --stats-file ${meta.id}.stats.csv \\
        --blocks-file ${meta.id}.blocks.tsv \\
        --summary-file ${meta.id}.summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$( hiphase -V | sed 's/hiphase //g')
    END_VERSIONS
    """
}
