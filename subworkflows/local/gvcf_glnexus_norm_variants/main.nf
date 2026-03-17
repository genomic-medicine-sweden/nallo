//
// Takes a channel of grouped gVCFs, merged them with GLNexus,
// adds FOUND_IN tag, normalizes and decomposes variants.
//

include { BCFTOOLS_PLUGINFIXPLOIDY                   } from '../../../modules/nf-core/bcftools/pluginfixploidy/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MULTISAMPLE } from '../../../modules/nf-core/bcftools/norm/main'
include { GLNEXUS                                    } from '../../../modules/nf-core/glnexus/main'
include { SENTIEON_GVCFTYPER                         } from '../../../modules/nf-core/sentieon/gvcftyper/main'
include { VCFEXPRESS                                  } from '../../../modules/nf-core/vcfexpress/main'

workflow GVCF_GLNEXUS_NORM_VARIANTS {
    take:
    ch_gvcfs                // channel: [mandatory] [ val(meta), path(gvcfs)     ]
    ch_tbis                 // channel: [mandatory] [ val(meta), path(tbis)      ]
    ch_bed                  // channel: [optional]  [ val(meta), path(input_bed) ]
    ch_fasta                // channel: [mandatory] [ val(meta), path(fasta)     ]
    ch_fai                  // channel: [mandatory] [ val(meta), path(fai)       ]
    variant_caller          // string: variant caller to tag the variants with, e.g. "deepvariant"
    ch_vcfexpress_prelude   // channel: [mandatory] [ val(meta), path(lua) ]

    main:
    ch_versions           = channel.empty()
    ch_merged_family_gvcf = channel.empty()

    if (variant_caller.equals("deepvariant")) {
        GLNEXUS(
            ch_gvcfs.map { meta, gvcfs -> [meta, gvcfs, []] },
            ch_bed,
        )

        ch_merged_family_gvcf = GLNEXUS.out.bcf
        ch_versions = ch_versions.mix(GLNEXUS.out.versions)

    } else if (variant_caller.equals("sentieon")) {

        ch_gvcfs
            .join(ch_tbis, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, gvcfs, tbis ->
                [meta, gvcfs, tbis, []]
            }
            .set { ch_gvcftyper_in }

        SENTIEON_GVCFTYPER(
            ch_gvcftyper_in,
            ch_fasta,
            ch_fai,
            [[], []],
            [[], []],
        )

        // Re-encoding haploid GTs to diploid  needs to occur _after_ GVCFtyper
        // in the gvcf joint-calling path, as the altered ploidy of haploid
        // to diploid GT no longer matches the sample PL field, which leads to a
        // crash in GVCFTyper
        BCFTOOLS_PLUGINFIXPLOIDY(
            SENTIEON_GVCFTYPER.out.vcf_gz.join(SENTIEON_GVCFTYPER.out.vcf_gz_tbi),
            [],
            [],
            [],
            []
        )

        ch_merged_family_gvcf = BCFTOOLS_PLUGINFIXPLOIDY.out.vcf

    }
    // Annotate with FOUND_IN tag - not sure what would happen if we do this before glnexus instead?
    ch_merged_family_gvcf
        .multiMap { meta, vcf ->
            vcf: [ meta, vcf ]
            sv_caller: meta.variant_caller
        }
        .set { ch_vcfexpress_input }

    ch_lua_file = ch_vcfexpress_prelude.map { meta, lua -> lua }

    VCFEXPRESS (
        ch_vcfexpress_input.vcf,
        ch_lua_file
    )

    // Decompose and normalize variants
    BCFTOOLS_NORM_MULTISAMPLE(
        VCFEXPRESS.out.vcf.map { meta, vcf -> [meta, vcf, []] },
        ch_fasta,
    )

    emit:
    vcf      = BCFTOOLS_NORM_MULTISAMPLE.out.vcf                                        // channel: [ val(meta), path(vcf) ]
    index    = BCFTOOLS_NORM_MULTISAMPLE.out.tbi.mix(BCFTOOLS_NORM_MULTISAMPLE.out.csi) // channel: [ val(meta), path(tbi/csi) ]
    versions = ch_versions                                                              // channel: [ path(versions.yml) ]
}
