include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'

workflow ALIGN_READS {

    take:
    ch_reads  // channel: [ val(meta), fastq ]
    ch_mmi    // channel: [ val(meta), mmi ]

    main:
    ch_versions = Channel.empty()

    MINIMAP2_ALIGN ( ch_reads, ch_mmi, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.csi)
        .set{ ch_bam_bai }

    emit:
    bam      = MINIMAP2_ALIGN.out.bam                              // channel: [ val(meta), bam ]
    bai      = MINIMAP2_ALIGN.out.csi                              // channel: [ val(meta), bai ]
    bam_bai  = MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.csi) // channel: [ val(meta), bam, bai ]
    versions = ch_versions                                         // channel: [ versions.yml ]
}

