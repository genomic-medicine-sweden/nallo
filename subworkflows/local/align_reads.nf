include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow ALIGN_READS {

    take:
    ch_fastq // channel: [ val(meta), fastq ]
    ch_mmi   // channel: [ val(meta), mmi ]

    main:
    ch_versions = Channel.empty()
    
    MINIMAP2_ALIGN( ch_fastq.combine( ch_mmi.map{ it[1] } ), true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    SAMTOOLS_INDEX ( MINIMAP2_ALIGN.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    MINIMAP2_ALIGN
        .out.bam
        .concat(SAMTOOLS_INDEX.out.bai)
        .groupTuple().flatten().collate(3)
        .set{ch_bam_bai}

    emit:
    bam_bai  = ch_bam_bai  // channel: [ val(meta), bam, bai]
    versions = ch_versions // channel: [ versions.yml ]
}

