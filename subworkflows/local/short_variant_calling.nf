// Is this the best way?
include { DEEPVARIANT } from '../../modules/local/google/deepvariant'
include { DEEPTRIO    } from '../../modules/local/google/deeptrio'
include { GLNEXUS     } from '../../modules/nf-core/glnexus'

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
        GLNEXUS ( DEEPVARIANT.out.gvcf.map { it [1] }.concat(ch_extra_gvcfs.map{ it[1] } ).collect().sort { it.name }.map{ [[id:"multisample"], it]} )
         
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
        ch_versions = ch_versions.mix(GLNEXUS.out.versions)
  
        ch_snp_calls_vcf = DEEPVARIANT.out.vcf
        ch_snp_calls_gvcf = DEEPVARIANT.out.gvcf
        ch_combined_bcf = GLNEXUS.out.bcf
    
    } else {
        ch_ped
            .combine(ch_bam_bai
                .map{ meta, bam, bai -> [meta.id, bam, bai]
                }
            )
            .branch{
                kid: it[0]['id']          == it[1]
                mom: it[0]['maternal_id'] == it[1]
                dad: it[0]['paternal_id'] == it[1]
            }
            .set{branch_result}

        branch_result
            .kid
            .join(branch_result.dad, remainder: true)
            .join(branch_result.mom, remainder: true)
            .map{[it[0], it[2], it[3], it[5], it[6], it[8], it[9]]} // [meta, kid_bam, kid_bai, dad_bam, dad_bai, mom_bam, mom_bai]
            .branch{ meta, kid_bam, kid_bai, dad_bam, dad_bai, mom_bam, mom_bai -> 
                is_trio: (dad_bam != null && dad_bai != null) && (mom_bam != null && mom_bai != null)
                    return tuple ( meta, kid_bam, kid_bai, dad_bam, dad_bai, mom_bam, mom_bai ) 
                no_trio: (dad_bam == null && dad_bai == null) || (mom_bam == null && mom_bai == null)
                    return tuple ( meta, kid_bam, kid_bai, [], [], [], [] )
            }
        .set{ch_samples}
        
        trio_kids = ch_samples.is_trio.map{[it[0], it[1], it[2]]}
        trio_dads = ch_samples.is_trio.map{[it[0], it[3], it[4]]}
        trio_moms = ch_samples.is_trio.map{[it[0], it[5], it[6]]}

        //non_trio_kids = ch_samples.no_trio.map{[it[0], it[1], it[2]]} // These can be someone elses mom or dad 
        //non_trio_dads = ch_samples.no_trio.map{[it[0], it[3], it[4]]} // These should all be empty
        //non_trio_moms = ch_samples.no_trio.map{[it[0], it[5], it[6]]} // These should all be empty
        
        DEEPTRIO( trio_kids.combine( ch_fasta.map{ it[1] } ).combine( ch_fai.map{ it[1] } ), trio_dads, trio_moms, 'PACBIO')
        
        // Get 3 vcf files back...we should collect all
        // Deal with multiple Parent VCFs in nextflow/GLNexus by naming them child.child/paternal/maternal 
        // and the person running multiple trios per family will have to deal with it later (me)
        // Do we even want to merge though?
        
        GLNEXUS ( DEEPTRIO.out.gvcf.map { it [1] }.collect().concat(ch_extra_gvcfs.map{ it[1] } ).collect().sort { it.name } )
    
        // Maaaybe do an is_parent check and if not, run DEEPVARIANT

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
