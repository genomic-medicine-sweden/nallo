process DNASCOPE_LONGREAD_CALL_SNVS {

   tag "$meta.id"
   label 'process_high'
   label 'sentieon'

   container "docker.io/clinicalgenomicslund/dnascope-longread:1.4.0-lrRPA-f9c8811"

   input:
   tuple val(meta),  path(bam), path(bai), path(interval_bed)
   tuple val(meta2), path(fasta)
   tuple val(meta3), path(fai)
   tuple val(meta4), path(model_bundle)
   tuple val(meta5), path(female_diploid_bed)
   tuple val(meta6), path(male_diploid_bed)
   tuple val(meta7), path(male_haploid_bed)
   val(tech)

   output:
   tuple val(meta), path("${prefix}.vcf.gz")      , emit: vcf
   tuple val(meta), path("${prefix}.vcf.gz.tbi")  , emit: vcf_tbi
   tuple val(meta), path("${prefix}.g.vcf.gz")    , emit: gvcf
   tuple val(meta), path("${prefix}.g.vcf.gz.tbi"), emit: gvcf_tbi
   path "versions.yml"                            , emit: versions

   script:
   prefix = task.ext.prefix ?: "${meta.id}"

   def diploid_intersect_cmd = ""
   def haploid_intersect_cmd = ""

   if(meta.sex == 1) {
       diploid_intersect_cmd = "bedtools intersect -a ${interval_bed} -b ${male_diploid_bed} > diploid_regions.bed"
       haploid_intersect_cmd = "bedtools intersect -a ${interval_bed} -b ${male_haploid_bed} > haploid_regions.bed"
   } else if(meta.sex == 2) {
       diploid_intersect_cmd = "bedtools intersect -a ${interval_bed} -b ${female_diploid_bed} > diploid_regions.bed"
   }

   """
   ${haploid_intersect_cmd}
   ${diploid_intersect_cmd}

   sentieon-cli dnascope-longread \\
        -t ${task.cpus} \\
        --tech ${tech} \\
        -r ${fasta} \\
        -i ${bam} \\
        -m ${model_bundle} \\
        ${diploid_intersect_cmd ? '--bed diploid_regions.bed' : "--bed ${interval_bed}"} \\
        ${haploid_intersect_cmd ? '--haploid_bed haploid_regions.bed' : ''} \\
        --gvcf \\
        --skip_mosdepth \\
        --skip_cnv \\
        --skip_svs \\
    ${prefix}.vcf.gz

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       sentieon-cli: 1.4.0-f9c8811
       sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
       bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
       samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
       bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
   END_VERSIONS
   """

   stub:
   prefix = task.ext.prefix ?: "${meta.id}"

   """
   touch ${prefix}.vcf.gz
   touch ${prefix}.vcf.gz.tbi
   touch ${prefix}.g.vcf.gz
   touch ${prefix}.g.vcf.gz.tbi

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       sentieon-cli: 1.4.0-f9c8811
       sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
       bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
       samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
       bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
   END_VERSIONS
   """

}
