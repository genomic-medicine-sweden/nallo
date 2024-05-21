include { MODKIT_PILEUP                                      } from '../../modules/local/modkit/pileup/main'
include { MODKIT_PILEUP as MODKIT_PILEUP_HAPLOTYPES          } from '../../modules/local/modkit/pileup/main'
include { TABIX_BGZIPTABIX as BGZIP_MODKIT_PILEUP            } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as BGZIP_MODKIT_PILEUP_HAPLOTYPES } from '../../modules/nf-core/tabix/bgziptabix/main'
workflow METHYLATION {

    take:
    ch_haplotagged_bam_bai // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_fai                 // channel: [ val(meta), fai ]
    ch_bed                 // channel: [ val(meta), bed ]

    main:
    ch_versions = Channel.empty()

    MODKIT_PILEUP(ch_haplotagged_bam_bai, ch_fasta, ch_fai, ch_bed)
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)

    MODKIT_PILEUP_HAPLOTYPES(ch_haplotagged_bam_bai, ch_fasta, ch_fai, ch_bed)
    ch_versions = ch_versions.mix(MODKIT_PILEUP_HAPLOTYPES.out.versions)

    BGZIP_MODKIT_PILEUP ( MODKIT_PILEUP.out.bed )

    MODKIT_PILEUP_HAPLOTYPES.out.haplotype_1
        .concat ( MODKIT_PILEUP_HAPLOTYPES.out.haplotype_2 )
        .concat ( MODKIT_PILEUP_HAPLOTYPES.out.ungrouped )
        .set { ch_bgzip_modkit_haplotypes_in }

    BGZIP_MODKIT_PILEUP_HAPLOTYPES ( ch_bgzip_modkit_haplotypes_in )

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

