//
// A subworkflow to plot binned zygosity and RHO-regions.
//

include { BCFTOOLS_ROH                              } from '../../../modules/nf-core/bcftools/roh/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_UNCOMPRESS } from '../../../modules/nf-core/bcftools/view/main'
include { RHOCALL_VIZ                               } from '../../../modules/nf-core/rhocall/viz/main'
include { CHROMOGRAPH as CHROMOGRAPH_AUTOZYG        } from '../../../modules/nf-core/chromograph/main'

workflow ANNOTATE_RHOCALLVIZ {
    take:
    ch_vcf // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_tbi // channel: [mandatory] [ val(meta), path(tbi) ]

    main:
    ch_versions = channel.empty()

    ch_vcf
        .join(ch_tbi, failOnMismatch:true, failOnDuplicate:true)
        .set { ch_vcf_tbi }

    BCFTOOLS_ROH(
        ch_vcf_tbi,
        [[], []],
        [],
        [],
        [],
        [],
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ROH.out.versions)

    BCFTOOLS_VIEW_UNCOMPRESS(
        ch_vcf_tbi,
        [],
        [],
        [],
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_UNCOMPRESS.out.versions)

    BCFTOOLS_VIEW_UNCOMPRESS.out.vcf
        .join(BCFTOOLS_ROH.out.roh)
        .multiMap { meta, vcf, roh ->
            vcf: [meta, vcf]
            roh: [meta, roh]
        }
        .set { ch_rhocall_viz_input }

    RHOCALL_VIZ(
        ch_rhocall_viz_input.vcf,
        ch_rhocall_viz_input.roh,
    )
    ch_versions = ch_versions.mix(RHOCALL_VIZ.out.versions)

    CHROMOGRAPH_AUTOZYG(
        RHOCALL_VIZ.out.bed,
        [[], []],
        [[], []],
        [[], []],
        [[], []],
        [[], []],
        [[], []],
    )
    ch_versions = ch_versions.mix(CHROMOGRAPH_AUTOZYG.out.versions)

    emit:
    versions          = ch_versions                   // channel: [ path(versions.yml) ]
    chromograph_plots = CHROMOGRAPH_AUTOZYG.out.plots // channel: [ val(meta), path(plot) ]
}
