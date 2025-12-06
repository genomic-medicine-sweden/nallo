//
// A subworkflow to plot binned coverage and zygosity regions.
//

include { BCFTOOLS_ROH                              } from '../../../modules/nf-core/bcftools/roh/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_UNCOMPRESS } from '../../../modules/nf-core/bcftools/view/main'
include { CHROMOGRAPH as RUN_CHROMOGRAPH            } from '../../../modules/nf-core/chromograph/main'
include { RHOCALL_VIZ                               } from '../../../modules/nf-core/rhocall/viz/main'
include { TIDDIT_COV                                } from '../../../modules/nf-core/tiddit/cov/main'

workflow CHROMOGRAPH {
    take:
    ch_bam_bai        // channel: [ val(meta), path(bam), path(bai) ]
    ch_vcf            // channel: [ val(meta), path(vcf) ]
    ch_tbi            // channel: [ val(meta), path(tbi) ]
    plot_coverage     // boolean
    plot_autozygosity // boolean

    main:
    ch_versions = channel.empty()
    ch_autozyg  = channel.of([[], []])
    ch_coverage = channel.of([[], []])


    if (plot_coverage) {
        TIDDIT_COV(
            ch_bam_bai,
            [[], []],
        )
        ch_versions = ch_versions.mix(TIDDIT_COV.out.versions)

        TIDDIT_COV.out.wig
            .map { meta, wig ->
                // To match the VCF meta, which only has ID due to being split by bcftools +split
                def new_meta = [id: meta.id]
                [new_meta, wig]
            }
            .set { ch_coverage }
    }

    if (plot_autozygosity) {
        ch_vcf
            .join(ch_tbi, failOnMismatch: true, failOnDuplicate: true)
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

        RHOCALL_VIZ.out.bed.set { ch_autozyg }
    }

    // Combine and filter only if there's data
    ch_autozyg
        .combine(ch_coverage)
        .filter { autozyg_meta, _autozyg, coverage_meta, _coverage ->
            if(!autozyg_meta || !coverage_meta)
                return true
            autozyg_meta.id == coverage_meta.id
        }
        .multiMap { autozyg_meta, autozyg, coverage_meta, coverage ->
            autozyg: [autozyg_meta, autozyg]
            coverage: [coverage_meta, coverage]
        }
        .set { ch_chromograph_input }

    RUN_CHROMOGRAPH(
        ch_chromograph_input.autozyg,
        ch_chromograph_input.coverage,
        [[], []],
        [[], []],
        [[], []],
        [[], []],
        [[], []],
    )
    ch_versions = ch_versions.mix(RUN_CHROMOGRAPH.out.versions)

    emit:
    chromograph_plots = RUN_CHROMOGRAPH.out.plots // channel: [ val(meta), path(plot) ]
    versions          = ch_versions               // channel: [ path(versions.yml) ]
}
