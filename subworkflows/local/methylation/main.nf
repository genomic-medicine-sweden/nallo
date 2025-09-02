include { MODKIT_PILEUP    } from '../../../modules/nf-core/modkit/pileup/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { PACBIO_METHYLATION_ANALYSIS } from '../pacbio_methylation_analysis/main'

workflow METHYLATION {

    take:
    ch_bam_bai             // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_bed                 // channel: [ val(meta), bed ]

    main:
    ch_versions = Channel.empty()

    // Performs pileups per haplotype if the phasing workflow is on, set in config
    MODKIT_PILEUP (ch_bam_bai, ch_fasta, ch_bed)
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)

    MODKIT_PILEUP.out.bed
        .transpose()
        .set { ch_bgzip_modkit_pileup_in }

    TABIX_BGZIPTABIX ( ch_bgzip_modkit_pileup_in )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    //
    // SUBWORKFLOW: Run PacBio methylation analysis with pb-cpg-tools and MethBat
    //
    PACBIO_METHYLATION_ANALYSIS (
        ch_bam_bai,
        ch_fasta,
        ch_bed
    )
    ch_versions = ch_versions.mix(PACBIO_METHYLATION_ANALYSIS.out.versions)

    emit:
    bed                 = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, bed, _tbi -> [ meta, bed ] } // channel: [ val(meta), path(bed) ]
    tbi                 = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, _bed, tbi -> [ meta, tbi ] } // channel: [ val(meta), path(tbi) ]
    cpg_scores          = PACBIO_METHYLATION_ANALYSIS.out.cpg_scores          // channel: [ val(meta), path(bed) ]
    cpg_pileup          = PACBIO_METHYLATION_ANALYSIS.out.cpg_pileup          // channel: [ val(meta), path(tsv) ]
    pmd_scores          = PACBIO_METHYLATION_ANALYSIS.out.pmd_scores          // channel: [ val(meta), path(txt) ]
    methylation_calls   = PACBIO_METHYLATION_ANALYSIS.out.methylation_calls   // channel: [ val(meta), path(tsv) ]
    methbat_summary     = PACBIO_METHYLATION_ANALYSIS.out.methbat_summary     // channel: [ val(meta), path(txt) ]
    versions            = ch_versions // channel: [ versions.yml ]
}
