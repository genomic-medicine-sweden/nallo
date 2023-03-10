/*
 * Structural variant calling
 */

include { DEEPVARIANT } from '../../modules/local/deepvariant.nf'
include { GLNEXUS } from '../../modules/local/glnexus.nf'

workflow SNP_CALLING {

  take:
  ch_bam_bai
  ch_fasta
  ch_fai
  ch_extra_gvcfs

  main:
  ch_snp_calls_vcf = Channel.empty()
  ch_snp_calls_gvcf = Channel.empty()
  ch_combined_bcf = Channel.empty()
  ch_extra_gvcfs = ch_extra_gvcfs.map{ it[0..1] }

  DEEPVARIANT ( ch_bam_bai.combine(ch_fasta).combine(ch_fai) )
  GLNEXUS ( DEEPVARIANT.out.snp_gvcf.collect().concat(ch_extra_gvcfs.map{it[1]}).collect() )
  
  ch_snp_calls_vcf = DEEPVARIANT.out.snp_vcf
  ch_snp_calls_gvcf = DEEPVARIANT.out.snp_gvcf
  ch_combined_bcf = GLNEXUS.out.snp_bcf

  emit:
  ch_snp_calls_vcf
  ch_snp_calls_gvcf
  ch_combined_bcf
}
