process HIFICNV {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::hificnv=0.1.7"
    container "quay.io/biocontainers/hificnv:0.1.7--h9ee0642_0"


    input:
    tuple val(meta), path(bam), path(bai), path(maf_vcf), path(expected_cn_bed)
    tuple val(meta2), path(fasta)
    path(exclude_bed)

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

    def expected_cn = expected_cn_bed ? "--expected-cn ${expected_cn_bed}" : ""
    def exclude = exclude_bed ? "--exclude ${exclude_bed}" : ""
    def maf = maf_vcf ? "--maf ${maf_vcf}" : ""

    if ("$maf_vcf" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
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
    mv hificnv.*.maf.bw ${prefix}.maf.bw
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
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.depth.bw
    touch ${prefix}.maf.bw
    touch ${prefix}.bedgraph
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$(echo \$(hificnv -V) | sed 's/hificnv //' )
    END_VERSIONS
    """
}

