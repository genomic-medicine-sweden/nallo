include { ADD_FOUND_IN_TAG                           } from '../../../modules/local/add_found_in_tag/main'
include { BCFTOOLS_CONCAT                             } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SINGLESAMPLE } from '../../../modules/nf-core/bcftools/norm/main'

//
// Workflow to concatenate and normalize variants
//
workflow VCF_CONCAT_NORM_VARIANTS {
    take:
    ch_vcfs        // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
    variant_caller // string: variant caller to tag the variants with, e.g. "deepvariant"

    main:
    ch_versions = Channel.empty()

    BCFTOOLS_CONCAT(
        ch_vcfs.map { meta, vcfs -> [meta, vcfs, []] },
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    // Annotate with FOUND_IN tag - not sure what would happen if we do this before glnexus instead?
    ADD_FOUND_IN_TAG(
        BCFTOOLS_CONCAT.out.vcf.map { meta, vcf -> [meta, vcf, []] },
        variant_caller,
    )
    ch_versions = ch_versions.mix(ADD_FOUND_IN_TAG.out.versions)

    BCFTOOLS_NORM_SINGLESAMPLE(
        ADD_FOUND_IN_TAG.out.vcf.map { meta, vcf -> [meta, vcf, []] },
        ch_fasta,
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM_SINGLESAMPLE.out.versions)

    emit:
    vcf                 = BCFTOOLS_NORM_SINGLESAMPLE.out.vcf                                         // channel: [ val(meta), path(vcf) ]
    index               = BCFTOOLS_NORM_SINGLESAMPLE.out.tbi.mix(BCFTOOLS_NORM_SINGLESAMPLE.out.csi) // channel: [ val(meta), path(tbi/csi) ]
    bcftools_concat_vcf = BCFTOOLS_CONCAT.out.vcf                                                    // channel: [ val(meta), path(vcf) ]
    versions            = ch_versions                                                                // channel: [ path(versions.yml) ]
}
