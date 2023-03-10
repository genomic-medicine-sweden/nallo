process DEEPVARIANT {
    tag "$sampleId"
    
    container 'google/deepvariant:1.5.0'
   
   publishDir 'data/interim/snp-calling/deepvariant', mode: 'copy'
    cpus = 16
    time '48h'

    input:
    tuple val(sampleId), path(reads), path(index), path(fasta), path(fai)

    output:
    path("${sampleId}.vcf.gz"),  emit: snp_vcf
    path("${sampleId}.vcf.gz"),  emit: snp_gvcf
    path("*.html"), emit: html

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --ref=${fasta} \\
        --reads=${reads} \\
        --model_type=PACBIO \\
        --output_vcf=${sampleId}.vcf.gz \\
        --output_gvcf=${sampleId}.gvcf.gz \\
        --num_shards=${task.cpus}
    """
}
