process HIPHASE {

    tag "$meta.id"
    label 'process_high'

    container "quay.io/biocontainers/hiphase:1.4.0--h9ee0642_0"

    input:
    tuple val(meta), path(vcfs), path(vcf_indices), path(bams), path(bais)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    val(output_bam)

    output:
    tuple val(meta), path("*.vcf.gz")      , emit: vcfs
    tuple val(meta), path("*.tbi")         , emit: vcfs_tbi     , optional: true
    tuple val(meta), path("*.csi")         , emit: vcfs_csi     , optional: true
    tuple val(meta), path("*.summary.tsv") , emit: summary_tsv  , optional: true
    tuple val(meta), path("*.summary.csv") , emit: summary_csv  , optional: true
    tuple val(meta), path("*.blocks.tsv")  , emit: blocks_tsv   , optional: true
    tuple val(meta), path("*.blocks.csv")  , emit: blocks_csv   , optional: true
    tuple val(meta), path("*.stats.tsv")   , emit: stats_tsv    , optional: true
    tuple val(meta), path("*.stats.csv")   , emit: stats_csv    , optional: true
    tuple val(meta), path("*.haplotag.tsv"), emit: haplotag_tsv , optional: true
    tuple val(meta), path("*.haplotag.csv"), emit: haplotag_csv , optional: true
    tuple val(meta), path("*.bam")         , emit: bams         , optional: true
    tuple val(meta), path("*.bam.bai")     , emit: bais         , optional: true
    tuple val(meta), path("*.bam.csi")     , emit: read_csis    , optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def vcfInputs = []
    def vcfOutputs = []
    def vcfNames = []
    for (vcf in vcfs) {
        vcfInputs.add('--vcf')
        vcfInputs.add(vcf)
        vcfOutputs.add('--output-vcf')
        vcfOutputs.add("${prefix}.vcf.gz")

        vcfNames.add(vcf.getName())
    }

    def bamInputs = []
    def bamOutputs = []
    def bamNames = []
    for (bam in bams) {
        bamInputs.add('--bam')
        bamInputs.add("${bam}")

        bamNames.add(bam.getName())

        if(output_bam) {
            bamOutputs.add('--output-bam')
            bamOutputs.add("${prefix}.bam")
        }
    }

    def uniqueVcfNames = new HashSet(vcfNames);
    if (uniqueVcfNames.size() < vcfNames.size()) {
        println("Name collision in input VCFs")
        exit 1
    }

    def uniqueBamNames = new HashSet(bamNames);
    if (uniqueBamNames.size() < bamNames.size()) {
        println("Name collision in input BAMs")
        exit 1
    }

    """
    hiphase \
        $args \
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        ${bamInputs.join(' ')} \\
        ${bamOutputs.join(' ')} \\
        ${vcfInputs.join(' ')} \\
        ${vcfOutputs.join(' ')}

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
