//
// Workflow to QC VCFs
//
include { DEEPVARIANT_VCFSTATSREPORT } from '../../../modules/nf-core/deepvariant/vcfstatsreport/main'
include { BCFTOOLS_STATS             } from '../../../modules/nf-core/bcftools/stats/main'

workflow QC_SNVS {
    take:
    ch_vcf              // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_normalized_vcf   // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_normalized_index // channel: [mandatory] [ val(meta), path(vcf) ]

    main:
    ch_versions = Channel.empty()

    DEEPVARIANT_VCFSTATSREPORT(
        ch_vcf
    )
    ch_versions = ch_versions.mix(DEEPVARIANT_VCFSTATSREPORT.out.versions)

    BCFTOOLS_STATS(
        ch_normalized_vcf.join(
            ch_normalized_index,
            failOnMismatch:true,
            failOnDuplicate:true
        ),
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    emit:
    vcfstatsreport = DEEPVARIANT_VCFSTATSREPORT.out.report // channel: [ val(meta), path(html) ]
    stats          = BCFTOOLS_STATS.out.stats              // channel: [ val(meta), path(txt) ]
    versions       = ch_versions                           // channel: [ path(versions.yml) ]
}
