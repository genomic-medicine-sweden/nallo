include { DEEPVARIANT                           } from '../../modules/local/google/deepvariant'
include { GLNEXUS                               } from '../../modules/nf-core/glnexus'
include { BCFTOOLS_VIEW_REGIONS                 } from '../../modules/local/bcftools/view_regions'
include { TABIX_TABIX as TABIX_EXTRA_GVCFS      } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DV               } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GLNEXUS          } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_DV } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_DV     } from '../../modules/nf-core/bcftools/sort/main'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_bam_bai
    ch_extra_gvcfs
    ch_fasta
    ch_fai
    ch_bed

    main:
    ch_snp_calls_vcf  = Channel.empty()
    ch_snp_calls_gvcf = Channel.empty()
    ch_combined_bcf   = Channel.empty()
    ch_versions       = Channel.empty()

    // Does splitting BAMs and copying to node make sense to reduce IO?

    // Only one of these is run depending on params.variant_caller (when clause condition is defined in the conf/modules.config)
    DEEPVARIANT               ( ch_bam_bai, ch_fasta, ch_fai )

    // Collect VCFs
    ch_snp_calls_vcf  = ch_snp_calls_vcf.mix(DEEPVARIANT.out.vcf)

    // Collect GVCFs
    ch_snp_calls_gvcf = ch_snp_calls_gvcf.mix(DEEPVARIANT.out.gvcf)

    // TODO: This only works with DeepVariant for now (remove PEPPER_MARGIN_DEEPVARIANT/Deeptrio?)

    // Extra gVCFs
    TABIX_EXTRA_GVCFS(ch_extra_gvcfs)

    ch_extra_gvcfs
        .join(TABIX_EXTRA_GVCFS.out.tbi)
        .groupTuple()
        .set{ ch_bcftools_view_regions_in }

    // This cuts all regions in BED file from extra gVCFS, better than nothing
    BCFTOOLS_VIEW_REGIONS( ch_bcftools_view_regions_in, ch_bed )

    // DV gVCFs
    TABIX_DV(ch_snp_calls_gvcf)

    ch_snp_calls_gvcf
        .groupTuple(size: params.parallel_snv)
        .join(TABIX_DV.out.tbi.groupTuple())
        .set{ bcftools_concat_dv_in }

    // Concat into one gVCF per sample & sort
    BCFTOOLS_CONCAT_DV ( bcftools_concat_dv_in )
    BCFTOOLS_SORT_DV   ( BCFTOOLS_CONCAT_DV.out.vcf )

    // Put DV and extra gvCFs together -> send to glnexus
    BCFTOOLS_SORT_DV.out.vcf
        .concat(BCFTOOLS_VIEW_REGIONS.out.vcf)
        .map { meta, gvcf -> [ ['id':'multisample'], gvcf ]}
        .groupTuple()
        .set{ ch_glnexus_in }

    // Multisample
    GLNEXUS( ch_glnexus_in, ch_bed )
    TABIX_GLNEXUS(GLNEXUS.out.bcf)

    // Get versions
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_REGIONS.out.versions)
    ch_versions = ch_versions.mix(TABIX_EXTRA_GVCFS.out.versions)
    ch_versions = ch_versions.mix(TABIX_DV.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_DV.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_DV.out.versions)
    ch_versions = ch_versions.mix(TABIX_GLNEXUS.out.versions)


    emit:
    snp_calls_vcf = BCFTOOLS_SORT_DV.out.vcf
    combined_bcf  = GLNEXUS.out.bcf
    versions      = ch_versions
}
