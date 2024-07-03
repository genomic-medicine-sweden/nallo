include { DEEPVARIANT                               } from '../../../modules/nf-core/deepvariant'
include { GLNEXUS                                   } from '../../../modules/nf-core/glnexus'
include { TABIX_TABIX as TABIX_DV                   } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DV_VCF               } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_DV     } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_DV_VCF } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_FILLTAGS                         } from '../../../modules/local/bcftools/filltags/main'
include { BCFTOOLS_NORM                             } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_DV         } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_DV_VCF     } from '../../../modules/nf-core/bcftools/sort/main'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_bam_bai_bed // channel: [mandatory] [ val(meta), path(bam), path(bai), path(call_region_bed) ]
    ch_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai         // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed         // channel: [optional] [ val(meta), path(input_bed) ]

    main:
    ch_snp_calls_vcf  = Channel.empty()
    ch_snp_calls_gvcf = Channel.empty()
    ch_combined_bcf   = Channel.empty()
    ch_versions       = Channel.empty()

    DEEPVARIANT ( ch_bam_bai_bed, ch_fasta, ch_fai, [[],[]] )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

    // Collect VCFs
    ch_snp_calls_vcf  = ch_snp_calls_vcf.mix(DEEPVARIANT.out.vcf)

    // Collect GVCFs
    ch_snp_calls_gvcf = ch_snp_calls_gvcf.mix(DEEPVARIANT.out.gvcf)

    // DV gVCFs
    TABIX_DV(ch_snp_calls_gvcf)
    ch_versions = ch_versions.mix(TABIX_DV.out.versions)

    ch_snp_calls_gvcf
        .groupTuple() // size not working here if there are less than specifed regions..
        .join(TABIX_DV.out.tbi.groupTuple())
        .set{ bcftools_concat_dv_in }

    // Concat into one gVCF per sample & sort
    BCFTOOLS_CONCAT_DV ( bcftools_concat_dv_in )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_DV.out.versions)

    BCFTOOLS_SORT_DV   ( BCFTOOLS_CONCAT_DV.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_DV.out.versions)

    // DV VCFs
    TABIX_DV_VCF(ch_snp_calls_vcf)
    ch_versions = ch_versions.mix(TABIX_DV_VCF.out.versions)

    ch_snp_calls_vcf
        .groupTuple() // size not working here if there are less than specifed regions..
        .join(TABIX_DV_VCF.out.tbi.groupTuple())
        .set{ bcftools_concat_dv_vcf_in }


    // Concat into one VCF per sample & sort
    BCFTOOLS_CONCAT_DV_VCF ( bcftools_concat_dv_vcf_in )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_DV_VCF.out.versions)

    BCFTOOLS_SORT_DV_VCF   ( BCFTOOLS_CONCAT_DV_VCF.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_DV_VCF.out.versions)

    // Put DV and extra gvCFs together -> send to glnexus
    BCFTOOLS_SORT_DV.out.vcf
        .map { meta, gvcf -> [ ['id':'multisample'], gvcf ]}
        .groupTuple()
        .set{ ch_glnexus_in }

    // Multisample
    GLNEXUS( ch_glnexus_in, ch_bed )
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)

    // Add allele count tag to multisample bcf
    BCFTOOLS_FILLTAGS ( GLNEXUS.out.bcf )
    ch_versions = ch_versions.mix(BCFTOOLS_FILLTAGS.out.versions)

    // Decompose and normalize variants
    BCFTOOLS_FILLTAGS.out.vcf
        .concat( BCFTOOLS_SORT_DV_VCF.out.vcf)
        .map { meta, vcf -> [ meta, vcf, [] ] }
        .set { bcftools_norm_in }

    BCFTOOLS_NORM ( bcftools_norm_in, ch_fasta )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    // Temporary solution while this workflow still outputs two types of vcfs
    BCFTOOLS_NORM.out.vcf
        .branch { meta, vcf ->
            multisample: meta.id == "multisample"
            singlesample: meta.id != "multisample"
        }
        .set { vcf_out }

    emit:
    snp_calls_vcf = vcf_out.singlesample // channel: [ val(meta), path(vcf) ]
    combined_bcf  = vcf_out.multisample  // channel: [ val(meta), path(bcf) ]
    versions      = ch_versions          // channel: [ path(versions.yml) ]
}
