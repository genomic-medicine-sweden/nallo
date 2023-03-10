/*
 * Structural variant calling
 */

include { SNIFFLES_SINGLE_SAMPLE } from '../../modules/local/sniffles_single_sample.nf'
include { SNIFFLES_MULTI_SAMPLE  } from '../../modules/local/sniffles_multi_sample.nf'

workflow STRUCTURAL_VARIANT_CALLING {

  take:
  ch_bam_bai
  ch_snfs

  main:
  ch_snfs = ch_snfs.map{ it[0..1] }
  ch_sv_calls_vcf = Channel.empty()

  SNIFFLES_SINGLE_SAMPLE( ch_bam_bai )
  SNIFFLES_MULTI_SAMPLE( SNIFFLES_SINGLE_SAMPLE.out.sv_snf.map{it[1]}.collect().concat(ch_snfs.map{it[1]}).collect())
  
  ch_sv_calls_vcf = SNIFFLES_MULTI_SAMPLE.out.multisample_vcf

  emit:
  ch_sv_calls_vcf
}
