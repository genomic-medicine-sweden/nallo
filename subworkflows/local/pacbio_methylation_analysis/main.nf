//
// PacBio DNA methylation analysis using pb-cpg-tools and MethBat
//

include { PBCPGTOOLS } from '../../../modules/local/pbcpgtools/main'
include { METHBAT    } from '../../../modules/local/methbat/main'

workflow PACBIO_METHYLATION_ANALYSIS {

    take:
    ch_bam_bai      // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference    // channel: [ val(meta), path(fasta) ]
    ch_regions      // channel: path(bed) - optional regions file

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run pb-cpg-tools for CpG score calculation
    //
    PBCPGTOOLS (
        ch_bam_bai,
        ch_reference,
        'aligned_bam_to_cpg_scores'
    )
    ch_versions = ch_versions.mix(PBCPGTOOLS.out.versions)

    //
    // MODULE: Run pb-cpg-tools for CpG pileup
    //
    PBCPGTOOLS_PILEUP = PBCPGTOOLS.cloneWithName('PBCPGTOOLS_PILEUP')
    PBCPGTOOLS_PILEUP (
        ch_bam_bai,
        ch_reference,
        'cpg_pileup'
    )

    //
    // MODULE: Run pb-cpg-tools for PMD calculation
    //
    PBCPGTOOLS_PMD = PBCPGTOOLS.cloneWithName('PBCPGTOOLS_PMD')
    PBCPGTOOLS_PMD (
        ch_bam_bai,
        ch_reference,
        'calculate_pmd'
    )

    //
    // MODULE: Run MethBat for methylation calling
    //
    METHBAT (
        ch_bam_bai,
        ch_reference,
        ch_regions
    )
    ch_versions = ch_versions.mix(METHBAT.out.versions)

    emit:
    cpg_scores          = PBCPGTOOLS.out.bed           // channel: [ val(meta), path(bed) ]
    cpg_pileup          = PBCPGTOOLS_PILEUP.out.tsv    // channel: [ val(meta), path(tsv) ]
    pmd_scores          = PBCPGTOOLS_PMD.out.txt       // channel: [ val(meta), path(txt) ]
    methylation_calls   = METHBAT.out.methylation_calls // channel: [ val(meta), path(tsv) ]
    methbat_summary     = METHBAT.out.summary          // channel: [ val(meta), path(txt) ]
    versions            = ch_versions                  // channel: [ versions.yml ]
}