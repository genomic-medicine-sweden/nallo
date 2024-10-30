include { MODKIT_PILEUP as MODKIT_PILEUP_UNPHASED } from '../../../modules/nf-core/modkit/pileup/main'
include { MODKIT_PILEUP as MODKIT_PILEUP_PHASED   } from '../../../modules/nf-core/modkit/pileup/main'
include { TABIX_BGZIPTABIX                        } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow METHYLATION {

    take:
    ch_bam_bai             // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_fai                 // channel: [ val(meta), fai ]
    ch_bed                 // channel: [ val(meta), bed ]
    phased                 // bool

    main:
    ch_versions = Channel.empty()

    if (phased) {
        MODKIT_PILEUP_PHASED (ch_bam_bai, ch_fasta, ch_bed)
        ch_versions = ch_versions.mix(MODKIT_PILEUP_PHASED.out.versions)

        MODKIT_PILEUP_PHASED.out.bed
            .transpose()
            .set { ch_bgzip_modkit_pileup_in }
    } else {
        MODKIT_PILEUP_UNPHASED (ch_bam_bai, ch_fasta, ch_bed)
        ch_versions = ch_versions.mix(MODKIT_PILEUP_UNPHASED.out.versions)

        ch_bgzip_modkit_pileup_in = MODKIT_PILEUP_UNPHASED.out.bed
    }

    TABIX_BGZIPTABIX ( ch_bgzip_modkit_pileup_in )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

