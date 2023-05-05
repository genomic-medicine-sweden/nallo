include { DEEPVARIANT } from '../../modules/local/deepvariant'
include { DEEPTRIO } from '../../modules/local/deeptrio'
include { GLNEXUS } from '../../modules/local/glnexus'
include { FIX_TRIO } from '../../modules/local/fix_trio'

workflow SHORT_VARIANT_CALLING {

  take:
  ch_bam_bai
  ch_extra_gvcfs
  ch_fasta
  ch_fai
  ch_ped

  main:
  ch_snp_calls_vcf = Channel.empty()
  ch_snp_calls_gvcf = Channel.empty()
  ch_combined_bcf = Channel.empty()
  ch_extra_gvcfs = Channel.empty()

  ch_versions = Channel.empty()

  if(!params.trio) {
    // First run DeepVariant with the aligned reads
    DEEPVARIANT ( ch_bam_bai.combine(ch_fasta.map { it[1] }).combine(ch_fai.map { it[1] }) )
  
    // Then run GlNexus to join-call genotypes (+ add previously run samples)
    GLNEXUS ( DEEPVARIANT.out.gvcf.map { it [1] }.concat(ch_extra_gvcfs.map{ it[1] } ).collect().sort { it.name } )
  
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)
  
    ch_snp_calls_vcf = DEEPVARIANT.out.vcf
    ch_snp_calls_gvcf = DEEPVARIANT.out.gvcf
    ch_combined_bcf = GLNEXUS.out.multisample_bcf
  
  } else {
    // Unfortuenate monstrosity: extract trios from PED file (TODO: add PED-validation)
    collected_bams = ch_bam_bai.map{ [it[1], it[2]] }.collect()
    // Combine the ped records with all of the bams..this might cause a lot of copying depending config..?
    combined = ch_ped.combine(collected_bams)
    // But unfortuenately i don't know how to do this in nextflow...
    FIX_TRIO ( ch_ped.combine(collected_bams.map{ [it] }) )
    // Get output back
    input_for_deeptrio = FIX_TRIO.out.trio.filter { it[0] =~ /0/ }.map{ [ it[1], it[2] ] }
    // Then run Deeptrio
    DEEPTRIO (input_for_deeptrio.combine(ch_fasta.map { it[1] }).combine(ch_fai.map { it[1] }).map{ [ it[0], it[2], it[3], it[1] ] } )
    // Get 3 vcf files back...we should collect all
    // Deal with multiple Parent VCFs in nextflow/GLNexus by naming them child.child/paternal/maternal 
    // and the person running multiple trios per family will have to deal with it later (me)
    GLNEXUS ( DEEPTRIO.out.gvcf.map { it [1] }.collect().concat(ch_extra_gvcfs.map{ it[1] } ).collect().sort { it.name } )
    
    ch_versions = ch_versions.mix(DEEPTRIO.out.versions)
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)
  
    ch_snp_calls_vcf = DEEPTRIO.out.vcf
    ch_snp_calls_gvcf = DEEPTRIO.out.gvcf
    ch_combined_bcf = GLNEXUS.out.multisample_bcf
  }
  


  emit:
  ch_snp_calls_vcf
  ch_snp_calls_gvcf
  ch_combined_bcf
  
  versions = ch_versions

}
