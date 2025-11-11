//
// Takes a channel of grouped gVCFs, merged them with GLNexus,
// adds FOUND_IN tag, normalizes and decomposes variants.
//
include { ADD_FOUND_IN_TAG                           } from '../../../modules/local/add_found_in_tag/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MULTISAMPLE } from '../../../modules/nf-core/bcftools/norm/main'
include { GLNEXUS                                    } from '../../../modules/nf-core/glnexus/main'
include { SENTIEON_GVCFTYPER                         } from '../../../modules/nf-core/sentieon/gvcftyper/main'

workflow GVCF_GLNEXUS_NORM_VARIANTS {
    take:
    ch_gvcfs       // channel: [mandatory] [ val(meta), path(gvcfs)     ]
    ch_gvcf_tbis   // channel: [mandatory] [ val(meta), path(tbis)      ]
    ch_bed         // channel: [optional]  [ val(meta), path(input_bed) ]
    ch_fasta       // channel: [mandatory] [ val(meta), path(fasta)     ]
    ch_fai         // channel: [mandatory] [ val(meta), path(fai)       ]
    variant_caller // string: variant caller to tag the variants with, e.g. "deepvariant"

    main:
    ch_versions           = Channel.empty()
    ch_merged_family_gvcf = Channel.empty()

    if(variant_caller.equals("deepvariant")) {
        GLNEXUS(
            ch_gvcfs.map { meta, gvcfs -> [meta, gvcfs, []] },
            ch_bed,
        )

        ch_merged_family_gvcf = GLNEXUS.out.bcf
        ch_versions = ch_versions.mix(GLNEXUS.out.versions)

    } else if(variant_caller.equals("sentieon")) {

        ch_gvcfs
            .join(ch_gvcf_tbis)
            .map{
                meta, gvcfs, tbis   ->
                [meta, gvcfs, tbis, []]
            }
            .set{ch_gvcftyper_in}

        SENTIEON_GVCFTYPER(
            ch_gvcftyper_in,
            ch_fasta,
            ch_fai,
            [[], []],
            [[], []]
        )

        ch_merged_family_gvcf = SENTIEON_GVCFTYPER.out.vcf_gz
        ch_versions = ch_versions.mix(SENTIEON_GVCFTYPER.out.versions)

    }
    // Annotate with FOUND_IN tag - not sure what would happen if we do this before glnexus instead?
    ADD_FOUND_IN_TAG(
        ch_merged_family_gvcf.map { meta, vcf -> [meta, vcf, []] },
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
