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
    ch_versions = ch_versions.mix(CRAMINO.out.versions)

    ch_bam_bai
        .combine( ch_bed.map { it[1] }.toList() ) // toList() enables passing [] if ch_bed is empty
        .set { mosdepth_in }

    MOSDEPTH ( mosdepth_in, ch_fasta )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

