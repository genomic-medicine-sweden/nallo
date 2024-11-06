//
// Annotate repeat expansions
//

include { BCFTOOLS_VIEW as COMPRESS_STRANGER } from '../../../modules/nf-core/bcftools/view/main'
include { STRANGER                           } from '../../../modules/nf-core/stranger/main'

workflow ANNOTATE_REPEAT_EXPANSIONS {
    take:
    ch_stranger_repeat_catalog // channel: [mandatory] [ path(stranger_repeat_catalog.json) ]
    ch_vcf             // channel: [mandatory] [ val(meta), path(vcf) ]

    main:
    ch_versions = Channel.empty()

    // Annotate, compress and index
    STRANGER ( ch_vcf, ch_stranger_repeat_catalog )
    ch_versions = ch_versions.mix(STRANGER.out.versions)

    COMPRESS_STRANGER (
        STRANGER.out.vcf.map { meta, vcf -> [meta, vcf, [] ] },
        [], [], []
    )
    ch_versions = ch_versions.mix(COMPRESS_STRANGER.out.versions)

    ch_vcf_idx = COMPRESS_STRANGER.out.vcf.join(COMPRESS_STRANGER.out.tbi, failOnMismatch:true, failOnDuplicate:true)

    emit:
    vcf_idx  = ch_vcf_idx   // channel: [ val(meta), path(vcf), path(tbi) ]
    versions = ch_versions  // channel: [ path(versions.yml) ]
}
