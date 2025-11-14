include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT   } from '../../../modules/nf-core/bcftools/sort/main'

//
// Workflow to concatenate and normalize variants
//
workflow VCF_CONCAT_SORT_VARIANTS {
    take:
    ch_vcfs_tbis // channel: [mandatory] [ val(meta), path(vcfs), path(tbis) ]

    main:
    ch_versions = channel.empty()

    BCFTOOLS_CONCAT(
        ch_vcfs_tbis
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    BCFTOOLS_SORT(
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    emit:
    vcf      = BCFTOOLS_SORT.out.vcf                            // channel: [ val(meta), path(vcf) ]
    index    = BCFTOOLS_SORT.out.tbi.mix(BCFTOOLS_SORT.out.csi) // channel: [ val(meta), path(tbi/csi) ]
    versions = ch_versions                                      // channel: [ path(versions.yml) ]
}
