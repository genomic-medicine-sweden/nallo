include { SAMTOOLS_FAIDX         } from '../../modules/nf-core/samtools/faidx/main'
include { MINIMAP2_INDEX         } from '../../modules/nf-core/minimap2/index/main'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/nf-core/gunzip/main'

workflow PREPARE_GENOME {

    take:
    fasta_in // channel: [ val(meta), fasta ]

    main:

    ch_versions = Channel.empty()

    // Will not catch cases where fasta is bgzipped
    if ( params.fasta.endsWith('.gz') ) {
        GUNZIP_FASTA(fasta_in)
            .gunzip
            .collect()
            .set{ch_fasta}

        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions.first())
    } else {
        fasta_in
            .set{ch_fasta}
    }

    SAMTOOLS_FAIDX ( ch_fasta )
    MINIMAP2_INDEX ( ch_fasta )

    // Gather versions
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())

    emit:
    fasta    = ch_fasta                           // channel: [ val(meta), fasta ]
    fai      = SAMTOOLS_FAIDX.out.fai.collect()   // channel: [ val(meta), fai ]
    mmi      = MINIMAP2_INDEX.out.index.collect() // channel: [ val(meta), mmi ]
    versions = ch_versions                        // channel: [ versions.yml ]
}
