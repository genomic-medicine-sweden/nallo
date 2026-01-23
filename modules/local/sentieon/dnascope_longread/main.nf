process DNASCOPE_LONGREAD_CALL_SNVS {

   tag "$meta.id"
   label 'process_high'
   label 'sentieon'

   container "docker.io/clinicalgenomicslund/dnascope-longread:1.4.0-lrRPA-f9c8811"

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
   tuple val("${task.process}"), val("sentieon-cli"), val("1.4.0-f9c8811"), emit: versions_sentieon_cli, topic: versions
   tuple val("${task.process}"), val("sentieon"), eval("sentieon driver --version 2>&1 | sed -e 's/sentieon-genomics-//g'"), emit: versions_sentieon, topic: versions
   tuple val("${task.process}"), val("bedtools"), eval("bedtools --version | sed -e 's/bedtools v//g'"), emit: versions_bedtools, topic: versions
   tuple val("${task.process}"), val("samtools"), eval("samtools --version 2>&1 | head -n 1 | sed 's/^.*samtools //; s/Using.*$//'"), emit: versions_samtools, topic: versions
   tuple val("${task.process}"), val("bcftools"), eval("bcftools --version 2>&1 | head -n 1 | sed 's/^.*bcftools //; s/ .*$//'"), emit: versions_bcftools, topic: versions

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
