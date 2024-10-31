include { MODKIT_PILEUP    } from '../../../modules/nf-core/modkit/pileup/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow METHYLATION {

    take:
    ch_bam_bai             // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_fai                 // channel: [ val(meta), fai ]
    ch_bed                 // channel: [ val(meta), bed ]

    main:
    ch_versions = Channel.empty()

    ch_bam_bai.dump()
    // Arguments to modkit are set according to whether phasing was performed in config
    MODKIT_PILEUP (ch_bam_bai, ch_fasta, ch_bed)
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)

    MODKIT_PILEUP.out.bed.dump()

    MODKIT_PILEUP.out.bed
        .transpose()
        .set { ch_bgzip_modkit_pileup_in }
    TABIX_BGZIPTABIX ( ch_bgzip_modkit_pileup_in )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    TABIX_BGZIPTABIX.out.dump()
    emit:
    bed      = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, bed, tbi -> [ meta, bed ] } // channel: [ val(meta), path(bed) ]
    tbi      = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, bed, tbi -> [ meta, tbi ] } // channel: [ val(meta), path(tbi) ]
    versions = ch_versions // channel: [ versions.yml ]
}

