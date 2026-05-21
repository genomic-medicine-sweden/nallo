//
// Takes a channel of grouped gVCFs, merged them with GLNexus,
// adds FOUND_IN tag, normalizes and decomposes variants.
//

include { BCFTOOLS_PLUGINFIXPLOIDY                   } from '../../../modules/nf-core/bcftools/pluginfixploidy/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MULTISAMPLE } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_MERGE                             } from '../../../modules/nf-core/bcftools/merge/main'
include { GLNEXUS                                    } from '../../../modules/nf-core/glnexus/main'
include { SENTIEON_GVCFTYPER                         } from '../../../modules/nf-core/sentieon/gvcftyper/main'
include { VCFEXPRESS                                 } from '../../../modules/nf-core/vcfexpress/main'

workflow GVCF_GLNEXUS_NORM_VARIANTS {
    take:
    ch_gvcfs                // channel: [mandatory] [ val(meta), path(gvcfs)     ]
    ch_tbis                 // channel: [mandatory] [ val(meta), path(tbis)      ]
    ch_bed                  // channel: [optional]  [ val(meta), path(input_bed) ]
    ch_fasta                // channel: [mandatory] [ val(meta), path(fasta)     ]
    ch_fai                  // channel: [mandatory] [ val(meta), path(fai)       ]
    ch_vcfexpress_prelude   // path: [mandatory] lua file

    main:
    ch_merged_family_gvcf = channel.empty()

    // Branching gVCFs and TBI channels by caller, as they need to be processed differently in the next steps
    ch_gvcfs
        .branch { meta, _gvcfs ->
            mitorsaw: meta.caller == "mitorsaw"
            sentieon: meta.caller == "sentieon"
            deepvariant: meta.caller == "deepvariant"
        }
        .set { branched_gvcfs }

    ch_tbis
        .branch { meta, _tbi ->
            mitorsaw: meta.caller == "mitorsaw"
            sentieon: meta.caller == "sentieon"
            deepvariant: meta.caller == "deepvariant"
        }
        .set { branched_tbis }

    // GLNEXUS processes deepvariant gVCFs
    GLNEXUS(
        branched_gvcfs.deepvariant.map { meta, gvcfs -> [meta, gvcfs, []] },
        ch_bed,
    )

    ch_merged_family_gvcf_glnexus = GLNEXUS.out.bcf

    // SENTIEON_GVCFTYPER processes sentieon gVCFs
    branched_gvcfs.sentieon
        .join(branched_tbis.sentieon, failOnMismatch: true, failOnDuplicate: true)
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

    // Re-encoding haploid GTs to diploid needs to occur _after_ GVCFtyper
    // in the gvcf joint-calling path, as the altered ploidy of haploid
    // to diploid GT no longer matches the sample PL field, which leads to a
    // crash in GVCFTyper
    BCFTOOLS_PLUGINFIXPLOIDY(
        SENTIEON_GVCFTYPER.out.vcf_gz.join(SENTIEON_GVCFTYPER.out.vcf_gz_tbi),
        [],
        [],
        [],
        [],
    )
    ch_merged_family_gvcf_sentieon = BCFTOOLS_PLUGINFIXPLOIDY.out.vcf

    // BCFTOOLS_MERGE is used for the mitochondrial specific caller
    branched_gvcfs.mitorsaw.join(branched_tbis.mitorsaw, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, gvcfs, tbis ->
            [meta, gvcfs, tbis, []]
        }
        .set { ch_bcftools_merge_in }

    ch_fasta.join(ch_fai).set { ch_fasta_fai }

    BCFTOOLS_MERGE(
        ch_bcftools_merge_in,
        ch_fasta_fai.collect(),
    )
    ch_merged_family_gvcf_bcftools = BCFTOOLS_MERGE.out.vcf

    ch_merged_family_gvcf_glnexus
        .mix(ch_merged_family_gvcf_sentieon)
        .mix(ch_merged_family_gvcf_bcftools)
        .set { ch_merged_family_gvcf }

    // Add FOUND_IN tag with VCFEXPRESS using the meta.caller information
    VCFEXPRESS(
        ch_merged_family_gvcf,
        ch_vcfexpress_prelude,
    )

    // Remove added caller information in meta
    VCFEXPRESS.out.vcf
        .map { meta, vcf ->
            [meta - meta.subMap('caller'), vcf, []]
        }
        .set { ch_bcftools_norm_input }


    // Decompose and normalize variants
    BCFTOOLS_NORM_MULTISAMPLE(
        ch_bcftools_norm_input,
        ch_fasta,
    )

    emit:
    vcf      = BCFTOOLS_NORM_MULTISAMPLE.out.vcf // channel: [ val(meta), path(vcf) ]
    index    = BCFTOOLS_NORM_MULTISAMPLE.out.tbi.mix(BCFTOOLS_NORM_MULTISAMPLE.out.csi) // channel: [ val(meta), path(tbi/csi) ]
}
