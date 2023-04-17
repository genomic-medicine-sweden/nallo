include { DEEPVARIANT } from '../../modules/local/deepvariant.nf'
include { GLNEXUS } from '../../modules/local/glnexus.nf'

workflow SHORT_VARIANT_CALLING {

  take:
  ch_bam_bai
  ch_extra_gvcfs
  ch_fasta
  ch_fai

  main:
  ch_snp_calls_vcf = Channel.empty()
  ch_snp_calls_gvcf = Channel.empty()
  ch_combined_bcf = Channel.empty()
  ch_extra_gvcfs = Channel.empty()

  ch_versions = Channel.empty()
  
  // First run DeepVariant with the aligned reads
  DEEPVARIANT ( ch_bam_bai.combine(ch_fasta.map { it[1] }).combine(ch_fai.map { it[1] }) )
  
  // Then run GlNexus to join-call genotypes (+ add previously run samples)
  GLNEXUS ( DEEPVARIANT.out.gvcf.map { it [1] }.concat(ch_extra_gvcfs.map{ it[1] } ).collect().sort { it.name } )
  
  ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
  ch_versions = ch_versions.mix(GLNEXUS.out.versions)

  ch_snp_calls_vcf = DEEPVARIANT.out.vcf
  ch_snp_calls_gvcf = DEEPVARIANT.out.gvcf
  ch_combined_bcf = GLNEXUS.out.multisample_bcf

  emit:
  ch_snp_calls_vcf
  ch_snp_calls_gvcf
  ch_combined_bcf
  
  versions = ch_versions

}
