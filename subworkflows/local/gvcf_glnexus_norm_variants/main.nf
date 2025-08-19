//
// Takes a channel of grouped gVCFs, merged them with GLNexus,
// adds FOUND_IN tag, normalizes and decomposes variants.
//
include { ADD_FOUND_IN_TAG                           } from '../../../modules/local/add_found_in_tag/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MULTISAMPLE } from '../../../modules/nf-core/bcftools/norm/main'
include { GLNEXUS                                    } from '../../../modules/nf-core/glnexus/main'

workflow GVCF_GLNEXUS_NORM_VARIANTS {
    take:
    ch_gvcfs       // channel: [mandatory] [ val(meta), path(gvcfs)]
    ch_bed         // channel: [optional]  [ val(meta), path(input_bed) ]
    ch_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
    variant_caller // string: variant caller to tag the variants with, e.g. "deepvariant"

    main:
    ch_versions = Channel.empty()

    GLNEXUS(
        ch_gvcfs,
        ch_bed,
    )
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)

    // Annotate with FOUND_IN tag - not sure what would happen if we do this before glnexus instead?
    ADD_FOUND_IN_TAG(
        GLNEXUS.out.bcf.map { meta, vcf -> [meta, vcf, []] },
        variant_caller,
    )
    ch_versions = ch_versions.mix(ADD_FOUND_IN_TAG.out.versions)

    // Decompose and normalize variants
    BCFTOOLS_NORM_MULTISAMPLE(
        ADD_FOUND_IN_TAG.out.vcf.map { meta, vcf -> [meta, vcf, []] },
        ch_fasta,
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM_MULTISAMPLE.out.versions)

    emit:
    vcf      = BCFTOOLS_NORM_MULTISAMPLE.out.vcf                                        // channel: [ val(meta), path(vcf) ]
    index    = BCFTOOLS_NORM_MULTISAMPLE.out.tbi.mix(BCFTOOLS_NORM_MULTISAMPLE.out.csi) // channel: [ val(meta), path(tbi/csi) ]
    versions = ch_versions                                                              // channel: [ path(versions.yml) ]
}
