include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MINIMAP2_ALIGN } from '../../modules/nf-core/samtools/index/main'

workflow ALIGN_READS {

    take:
    ch_fastq // channel: [ val(meta), fastq ]
    ch_mmi   // channel: [ val(meta), mmi ]

    main:
    ch_versions = Channel.empty()

    // Remap index
    ch_mmi = ch_mmi.map{ it [1] }

    
    MINIMAP2_ALIGN( ch_fastq.combine( ch_mmi ), true, false, false )
    SAMTOOLS_INDEX_MINIMAP2_ALIGN ( MINIMAP2_ALIGN.out.bam )
    
    MINIMAP2_ALIGN.out.bam
        .join(SAMTOOLS_INDEX_MINIMAP2_ALIGN.out.bai)
        .set{ ch_bam_bai }

    // Gather versions
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MINIMAP2_ALIGN.out.versions.first())

    emit:
    bam      = MINIMAP2_ALIGN.out.bam                // channel: [ val(meta), bam ]
    bai      = SAMTOOLS_INDEX_MINIMAP2_ALIGN.out.bai // channel: [ val(meta), bai ]
    bam_bai  = ch_bam_bai                            // channel: [ val(meta), bam, bai]
    versions = ch_versions                           // channel: [ versions.yml ]
}

