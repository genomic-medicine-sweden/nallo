include { MODKIT_PILEUP    } from '../../../modules/nf-core/modkit/pileup/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow METHYLATION {

    take:
    ch_bam_bai             // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_fai                 // channel: [ val(meta), fai ]
    ch_bed                 // channel: [ val(meta), bed ]
    phased                 // bool

    main:
    ch_versions = Channel.empty()

    // Arguments to modkit are set according to whether phasing was performed in config
    MODKIT_PILEUP (ch_bam_bai, ch_fasta, ch_bed)
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)

    // Output is different between phased and unphased modes.
    // We must deal with that here to get consistent subworkflow output
    if (phased) {
        MODKIT_PILEUP.out.bed
            .transpose()
            .set { ch_bgzip_modkit_pileup_in }
    } else {
        ch_bgzip_modkit_pileup_in = MODKIT_PILEUP.out.bed
    }

    TABIX_BGZIPTABIX ( ch_bgzip_modkit_pileup_in )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

