process DNASCOPE_LONGREAD_CALL_SNVS {

   tag "$meta.id"
   label 'process_high'
   label 'sentieon'

   container "docker.io/clinicalgenomicslund/dnascope-longread:1.5.1"

   input:
   tuple val(meta),  path(bam), path(bai), path(diploid_intervals_bed), path(haploid_intervals_bed)
   tuple val(meta2), path(fasta)
   tuple val(meta3), path(fai)
   tuple val(meta4), path(model_bundle)
   val(tech)

   output:
   tuple val(meta), path("${prefix}.vcf.gz")      , emit: vcf
   tuple val(meta), path("${prefix}.vcf.gz.tbi")  , emit: vcf_tbi
   tuple val(meta), path("${prefix}.g.vcf.gz")    , emit: gvcf
   tuple val(meta), path("${prefix}.g.vcf.gz.tbi"), emit: gvcf_tbi
   tuple val("${task.process}"), val("sentieon-cli"), eval("echo '1.4.0-f9c8811'"), topic: versions, emit: versions_sentieon_cli
   tuple val("${task.process}"), val("sentieon"), eval("sentieon driver --version 2>&1 | sed -e 's/sentieon-genomics-//g'"), topic: versions, emit: versions_sentieon
   tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools
   tuple val("${task.process}"), val("samtools"), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
   tuple val("${task.process}"), val("bcftools"), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

   script:
   prefix = task.ext.prefix ?: "${meta.id}"

   def haploid_bed_arg = haploid_intervals_bed ? "--haploid_bed ${haploid_intervals_bed}" : ""

   """
   sentieon-cli dnascope-longread \\
        -t ${task.cpus} \\
        --tech ${tech} \\
        -r ${fasta} \\
        -i ${bam} \\
        -m ${model_bundle} \\
        --bed ${diploid_intervals_bed} \\
        ${haploid_bed_arg} \\
        --gvcf \\
        --skip_mosdepth \\
        --skip_cnv \\
        --skip_svs \\
    ${prefix}.vcf.gz

   """

   stub:
   prefix = task.ext.prefix ?: "${meta.id}"
   """
   echo "" | bgzip > ${prefix}.vcf.gz
   touch ${prefix}.vcf.gz.tbi
   echo "" | bgzip > ${prefix}.g.vcf.gz
   touch ${prefix}.g.vcf.gz.tbi
   """

}
