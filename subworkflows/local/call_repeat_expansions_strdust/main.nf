include { STRDUST          } from '../../../modules/nf-core/strdust/'
include { ADD_FOUND_IN_TAG } from '../../../modules/local/add_found_in_tag/main'
include { BCFTOOLS_MERGE   } from '../../../modules/nf-core/bcftools/merge/'

workflow CALL_REPEAT_EXPANSIONS_STRDUST {
    take:
    ch_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta   // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai     // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed     // channel: [mandatory] [ val(meta), path(bed) ]

    main:
    ch_versions = Channel.empty()

    STRDUST(
        ch_bam_bai,
        ch_fasta,
        ch_fai,
        ch_bed,
    )
    ch_versions.mix(STRDUST.out.versions)

    ADD_FOUND_IN_TAG(
        STRDUST.out.vcf.join(STRDUST.out.tbi),
        "STRdust",
    )

    ADD_FOUND_IN_TAG.out.vcf
        .join(ADD_FOUND_IN_TAG.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi -> [[id: meta.family_id], vcf, tbi] }
        .groupTuple()
        .set { ch_bcftools_merge_in }

    BCFTOOLS_MERGE(
        ch_bcftools_merge_in,
        [[], []],
        [[], []],
        [[], []],
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    emit:
    sample_vcf = STRDUST.out.vcf // channel: [ val(meta), path(vcf) ]
    sample_tbi = STRDUST.out.tbi // channel: [ val(meta), path(tbi) ]
    family_vcf = BCFTOOLS_MERGE.out.vcf // channel: [ val(meta), path(vcf) ]
    family_tbi = BCFTOOLS_MERGE.out.index // channel: [ val(meta), path(tbi) ]
    versions   = ch_versions // channel: [ versions.yml ]
}
