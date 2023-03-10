process GLNEXUS {
    tag "$sampleId"
    container 'quay.io/mlin/glnexus:v1.2.7'
    publishDir 'data/interim/snp-calling/glnexus', mode: 'copy'
    cpus = 16
    time '48h'

    input:
    path(gvcfs)

    output:
    path("*.vcf.gz"),  emit: snp_bcf

    """
    glnexus_cli \\
      --config DeepVariant_unfiltered \\
      ${gvcfs} \\
      --threads ${task.cpus} \\
      | bcftools view - \\
      | bgzip -@ ${task.cpus} -c \\
      > multisample.vcf.gz
    """
}
