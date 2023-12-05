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
    ch_bam_bai
        .combine(ch_bed.map{ meta, bed -> bed })
        //.map{ meta, bam, bai -> [ meta, bam, bai, [] ] }
        .set{ ch_mosdepth_in }

    CRAMINO (ch_bam_bai)
    MOSDEPTH(ch_mosdepth_in, ch_fasta)

    // Gather versions
    ch_versions = ch_versions.mix(CRAMINO.out.versions.first())
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

