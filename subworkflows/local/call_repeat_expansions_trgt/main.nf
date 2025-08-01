include { TRGT_GENOTYPE    } from '../../../modules/nf-core/trgt/genotype/main'
include { SAMTOOLS_INDEX   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT    } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'
include { ADD_FOUND_IN_TAG } from '../../../modules/local/add_found_in_tag/main'
include { BCFTOOLS_SORT    } from '../../../modules/nf-core/bcftools/sort/main'
include { TRGT_MERGE       } from '../../../modules/nf-core/trgt/merge/main'
include { BCFTOOLS_INDEX   } from '../../../modules/nf-core/bcftools/index/main'

workflow CALL_REPEAT_EXPANSIONS_TRGT {
    take:
    ch_bam_bai  // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai      // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed      // channel: [mandatory] [ val(meta), path(bed) ]
    cram_output // bool: Publish alignments as CRAM (true) or BAM (false)

    main:

    ch_versions = Channel.empty()

    ch_bam_bai
        .map { meta, bam, bai -> [meta, bam, bai, meta.sex == 1 ? 'XY' : 'XX'] }
        .set { ch_trgt_input }

    // Run TRGT
    TRGT_GENOTYPE(
        ch_trgt_input,
        ch_fasta,
        ch_fai,
        ch_bed,
    )
    ch_versions = ch_versions.mix(TRGT_GENOTYPE.out.versions)

    // Sort and index bam
    SAMTOOLS_SORT(
        TRGT_GENOTYPE.out.bam,
        [[], []],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Publish spanning reads as CRAM if requested
    if (cram_output) {
        SAMTOOLS_CONVERT(
            SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true),
            ch_fasta,
            ch_fai,
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)
    }

    // Add FOUND_IN=TRGT tag
    ADD_FOUND_IN_TAG(
        TRGT_GENOTYPE.out.vcf.map { meta, vcf -> [meta, vcf, []] },
        "TRGT",
    )

    // Sort and index bcf
    BCFTOOLS_SORT(ADD_FOUND_IN_TAG.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    BCFTOOLS_SORT.out.vcf
        .join(BCFTOOLS_SORT.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, bcf, csi -> [[id: meta.family_id], bcf, csi] }
        .groupTuple()
        .set { ch_trgt_merge_in }

    TRGT_MERGE(
        ch_trgt_merge_in,
        [[], []],
        [[], []],
    )
    ch_versions = ch_versions.mix(TRGT_MERGE.out.versions)

    BCFTOOLS_INDEX(
        TRGT_MERGE.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

    emit:
    sample_vcf = BCFTOOLS_SORT.out.vcf // channel: [ val(meta), path(vcf) ]
    sample_tbi = BCFTOOLS_SORT.out.tbi // channel: [ val(meta), path(tbi) ]
    family_vcf = TRGT_MERGE.out.vcf // channel: [ val(meta), path(vcf) ]
    family_tbi = BCFTOOLS_INDEX.out.tbi // channel: [ val(meta), path(tbi) ]
    sample_bam = SAMTOOLS_SORT.out.bam // channel: [ val(meta), path(bam) ]
    sample_bai = SAMTOOLS_INDEX.out.bai // channel: [ val(meta), path(bai) ]
    versions   = ch_versions // channel: [ versions.yml ]
}
