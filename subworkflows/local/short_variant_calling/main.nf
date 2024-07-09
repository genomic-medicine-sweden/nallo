include { DEEPVARIANT                             } from '../../../modules/nf-core/deepvariant'
include { GLNEXUS                                 } from '../../../modules/nf-core/glnexus'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_GVCF } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF  } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_FILLTAGS                       } from '../../../modules/local/bcftools/filltags/main'
include { BCFTOOLS_NORM                           } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_GVCF     } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_VCF      } from '../../../modules/nf-core/bcftools/sort/main'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_bam_bai_bed // channel: [mandatory] [ val(meta), path(bam), path(bai), path(call_region_bed) ]
    ch_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai         // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed         // channel: [optional] [ val(meta), path(input_bed) ]

    main:
    ch_versions = Channel.empty()

    DEEPVARIANT ( ch_bam_bai_bed, ch_fasta, ch_fai, [[],[]] )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

    // gVCF
    DEEPVARIANT.out.gvcf
        .map { meta, vcf -> [ groupKey(meta, meta.num_intervals ), vcf ] }
        .groupTuple()
        .join( DEEPVARIANT.out.gvcf_tbi
            .map { meta, vcf -> [ groupKey(meta, meta.num_intervals ), vcf ] }
            .groupTuple()
        )
        .set{ bcftools_concat_gvcf_in }

    // Concat into one gVCF per sample & sort
    BCFTOOLS_CONCAT_GVCF ( bcftools_concat_gvcf_in )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_GVCF.out.versions)

    BCFTOOLS_SORT_GVCF ( BCFTOOLS_CONCAT_GVCF.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_GVCF.out.versions)

    // VCF
    DEEPVARIANT.out.vcf
        .map { meta, vcf -> [ groupKey(meta, meta.num_intervals ), vcf ] }
        .groupTuple()
        .join( DEEPVARIANT.out.vcf_tbi
            .map { meta, vcf -> [ groupKey(meta, meta.num_intervals ), vcf ] }
            .groupTuple()
        )
        .set{ bcftools_concat_vcf_in }

    // Concat into one VCF per sample & sort
    BCFTOOLS_CONCAT_VCF ( bcftools_concat_vcf_in )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_VCF.out.versions)

    BCFTOOLS_SORT_VCF ( BCFTOOLS_CONCAT_VCF.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_VCF.out.versions)

    // Multisample
    BCFTOOLS_SORT_GVCF.out.vcf
        .map { meta, gvcf -> [ ['id':'multisample'], gvcf ] }
        .groupTuple()
        .set{ glnexus_in }

    GLNEXUS( glnexus_in, ch_bed )
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)

    // Add allele count tag to multisample bcf
    BCFTOOLS_FILLTAGS ( GLNEXUS.out.bcf )
    ch_versions = ch_versions.mix(BCFTOOLS_FILLTAGS.out.versions)

    // Decompose and normalize variants
    BCFTOOLS_FILLTAGS.out.vcf
        .concat( BCFTOOLS_SORT_VCF.out.vcf)
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
