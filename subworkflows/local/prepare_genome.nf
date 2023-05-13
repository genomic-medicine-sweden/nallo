include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'

workflow PREPARE_GENOME {

    take:
    ch_fasta // channel: [ val(meta), fasta ]
    
    main:
    ch_versions = Channel.empty()

    SAMTOOLS_FAIDX ( ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    MINIMAP2_INDEX ( ch_fasta )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())
    
    emit:
    fasta = ch_fasta                  // channel: [ val(meta), fasta ]
    fai   = SAMTOOLS_FAIDX.out.fai    // channel: [ val(meta), fai ]
    mmi   = MINIMAP2_INDEX.out.index  // channel: [ val(meta), mmi ]

    versions = ch_versions            // channel: [ versions.yml ]
}
