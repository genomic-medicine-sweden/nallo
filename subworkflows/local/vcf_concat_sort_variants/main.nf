include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT   } from '../../../modules/nf-core/bcftools/sort/main'

//
// Workflow to concatenate and sort variants
//
workflow VCF_CONCAT_SORT_VARIANTS {
    take:
    ch_vcfs_tbis // channel: [mandatory] [ val(meta), path(vcfs), path(tbis) ]

    main:
    BCFTOOLS_CONCAT(
        ch_vcfs_tbis
    )

    BCFTOOLS_SORT(
        BCFTOOLS_CONCAT.out.vcf
    )

    emit:
    vcf      = BCFTOOLS_SORT.out.vcf                            // channel: [ val(meta), path(vcf) ]
    index    = BCFTOOLS_SORT.out.tbi.mix(BCFTOOLS_SORT.out.csi) // channel: [ val(meta), path(tbi/csi) ]
}
