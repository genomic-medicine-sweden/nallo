include { MODKIT_PILEUP                             } from '../../modules/local/modkit/pileup/main'
include { MODKIT_PILEUP as MODKIT_PILEUP_HAPLOTYPES } from '../../modules/local/modkit/pileup/main'

workflow METHYLATION {

    take:
    ch_haplotagged_bam_bai // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_fai                 // channel: [ val(meta), fai ]
    ch_bed                 // channel: [ val(meta), bed ]

    main:
    ch_versions = Channel.empty()

    MODKIT_PILEUP(ch_haplotagged_bam_bai, ch_fasta, ch_fai, ch_bed)
    MODKIT_PILEUP_HAPLOTYPES(ch_haplotagged_bam_bai, ch_fasta, ch_fai, ch_bed)

    // Get versions
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)
    ch_versions = ch_versions.mix(MODKIT_PILEUP_HAPLOTYPES.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

