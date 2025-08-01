process HIPHASE {

    tag "${meta.id}"
    label 'process_high'

    container "quay.io/biocontainers/hiphase:1.4.0--h9ee0642_0"

    input:
    tuple val(meta), path(vcfs), path(vcf_indices), path(bams), path(bais)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    val output_bam

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcfs
    tuple val(meta), path("*.tbi"), emit: vcfs_tbi, optional: true
    tuple val(meta), path("*.csi"), emit: vcfs_csi, optional: true
    tuple val(meta), path("*.summary.tsv"), emit: summary_tsv, optional: true
    tuple val(meta), path("*.summary.csv"), emit: summary_csv, optional: true
    tuple val(meta), path("*.blocks.tsv"), emit: blocks_tsv, optional: true
    tuple val(meta), path("*.blocks.csv"), emit: blocks_csv, optional: true
    tuple val(meta), path("*.stats.tsv"), emit: stats_tsv, optional: true
    tuple val(meta), path("*.stats.csv"), emit: stats_csv, optional: true
    tuple val(meta), path("*.haplotag.tsv"), emit: haplotag_tsv, optional: true
    tuple val(meta), path("*.haplotag.csv"), emit: haplotag_csv, optional: true
    tuple val(meta), path("*.bam"), emit: bams, optional: true
    tuple val(meta), path("*.bam.bai"), emit: bais, optional: true
    tuple val(meta), path("*.bam.csi"), emit: read_csis, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def bamNames = []
    def vcfNames = []

    def vcf_args = vcfs
        .collectMany { file ->
            [
                "--vcf",
                file,
                "--output-vcf",
                "${prefix}_phased.vcf.gz",
            ]
        }
        .join(" ")

    def bam_args = bams
        .collectMany { file ->
            [
                "--bam",
                file,
                output_bam ? '--output-bam' : '',
                output_bam ? "${prefix}_haplotagged.bam" : '',
            ]
        }
        .join(" ")

    vcfs.each { vcf ->
        vcfNames.add(vcf.getName())
    }

    bams.each { bam ->
        bamNames.add(bam.getName())
    }

    def uniqueVcfNames = new Set(vcfNames)
    if (uniqueVcfNames.size() < vcfNames.size()) {
        println("Name collision in input VCFs")
        exit(1)
    }

    def uniqueBamNames = new Set(bamNames)
    if (uniqueBamNames.size() < bamNames.size()) {
        println("Name collision in input BAMs")
        exit(1)
    }

    """
    hiphase \
        ${args} \
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        ${bam_args} \\
        ${vcf_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$( hiphase -V | sed 's/hiphase //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$( hiphase -V | sed 's/hiphase //g')
    END_VERSIONS
    """
}
