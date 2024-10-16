process HIFICNV {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::hificnv=1.0.0"
    container "quay.io/biocontainers/hificnv:1.0.0--h9ee0642_0"

    input:
    tuple val(meta), path(bam), path(bai), path(maf_vcf), val(sex)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(expected_xy_bed)
    tuple val(meta4), path(expected_xx_bed)
    tuple val(meta5), path(exclude_bed)

    output:
    tuple val(meta), path("*.vcf.gz")  , emit: vcf
    tuple val(meta), path("*.depth.bw"), emit: depth_bw
    tuple val(meta), path("*.maf.bw")  , emit: maf_bw, optional: true
    tuple val(meta), path("*.bedgraph"), emit: cnval
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    def expected_cn = sex == 1 ? "--expected-cn ${expected_xy_bed}" : sex == 2 ? "--expected-cn ${expected_xx_bed}" : ""
    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ""
    def maf = maf_vcf ? "--maf ${maf_vcf}" : ""

    if ("$maf_vcf" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    mv_maf = maf ? "mv hificnv.*.maf.bw ${prefix}.maf.bw" : ''

    """
    hificnv \\
        $args \\
        --bam ${bam} \\
        $expected_cn \\
        $exclude \\
        $maf \\
        --ref ${fasta} \\
        --threads ${task.cpus}

    mv hificnv.*.vcf.gz ${prefix}.vcf.gz
    mv hificnv.*.depth.bw ${prefix}.depth.bw
    $mv_maf
    mv hificnv.*.copynum.bedgraph ${prefix}.copynum.bedgraph
    mv *.log ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(echo \$(hificnv -V) | sed 's/hificnv //' )
    END_VERSIONS
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"

    if ("$maf_vcf" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    def maf = maf_vcf ? "touch ${prefix}.maf.bw" : ""
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.depth.bw
    $maf
    touch ${prefix}.bedgraph
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(echo \$(hificnv -V) | sed 's/hificnv //' )
    END_VERSIONS
    """
}

