include { CRAMINO  } from '../../modules/local/cramino'
include { MOSDEPTH } from '../../modules/nf-core/mosdepth/main.nf'

workflow QC_ALIGNED_READS {

    take:
    ch_bam_bai // channel: [ val(meta), [bam, bai] ]
    ch_fasta   // channel: [ val(meta), fasta ]
    ch_bed     // channel: [ val(meta), bed ]

    main:
    ch_versions = Channel.empty()

    // Prepare inputs

    CRAMINO (ch_bam_bai)
    MOSDEPTH(ch_bam_bai, ch_bed.map{ it[1] }, ch_fasta)

    // Gather versions
    ch_versions = ch_versions.mix(CRAMINO.out.versions.first())
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

