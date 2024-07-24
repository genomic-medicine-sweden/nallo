include { MODKIT_PILEUP as MODKIT_PILEUP_UNPHASED          } from '../../modules/nf-core/modkit/pileup/main'
include { MODKIT_PILEUP as MODKIT_PILEUP_PHASED            } from '../../modules/nf-core/modkit/pileup/main'
include { TABIX_BGZIPTABIX as BGZIP_MODKIT_PILEUP_UNPHASED } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as BGZIP_MODKIT_PILEUP_PHASED   } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow METHYLATION {

    take:
    ch_haplotagged_bam_bai // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_fai                 // channel: [ val(meta), fai ]
    ch_bed                 // channel: [ val(meta), bed ]

    main:
    ch_versions = Channel.empty()

    // Run modkit pileup once without dividing by HP-tag and once with
    MODKIT_PILEUP_UNPHASED (ch_haplotagged_bam_bai, ch_fasta, ch_bed)
    ch_versions = ch_versions.mix(MODKIT_PILEUP_UNPHASED.out.versions)

    MODKIT_PILEUP_PHASED (ch_haplotagged_bam_bai, ch_fasta, ch_bed)
    ch_versions = ch_versions.mix(MODKIT_PILEUP_PHASED.out.versions)

    // Bgzip and index output "BED"
    BGZIP_MODKIT_PILEUP_UNPHASED ( MODKIT_PILEUP_UNPHASED.out.bed )
    ch_versions = ch_versions.mix(BGZIP_MODKIT_PILEUP_UNPHASED.out.versions)

    MODKIT_PILEUP_PHASED.out.bed
        .transpose()
        .set { ch_bgzip_modkit_pileup_phased_in }

    BGZIP_MODKIT_PILEUP_PHASED ( ch_bgzip_modkit_pileup_phased_in )
    ch_versions = ch_versions.mix(BGZIP_MODKIT_PILEUP_PHASED.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

