include { MODKIT_PILEUP                             } from '../../modules/local/modkit/pileup/main'
include { MODKIT_PILEUP as MODKIT_PILEUP_HAPLOTYPES } from '../../modules/local/modkit/pileup/main'

workflow METHYLATION {

    take:
    ch_haplotagged_bam_bai // channel: [ val(meta), bam, bai ]
    ch_vcf     // channel: [ val(meta), vcf.gz ]
    ch_fasta   // channel: [ val(meta), fasta ]
    ch_fai     // channel: [ val(meta), fai ]

    main:
    ch_versions = Channel.empty()
    
    
    // Sometimes you may not want per haplotype methylation (e.g. chrX males)
    // ...Maybe we could to this only on male chrX and chrY?
    // TODO: Fix inconsistent naming between hp1/2 and haplotagged bams?
    
    MODKIT_PILEUP( ch_haplotagged_bam_bai, ch_fasta, ch_fai)
    MODKIT_PILEUP_HAPLOTYPES( ch_haplotagged_bam_bai, ch_fasta, ch_fai)
    
    // Sometimes you may want per haplotype methylation
    
    // Get versions
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)
    ch_versions = ch_versions.mix(MODKIT_PILEUP_HAPLOTYPES.out.versions)
    
    emit:
    //ch_combined_bedmethyl = MODKIT_PILEUP.out.bed // channel: [ val(meta), bed.gz, tbi ]
    versions              = ch_versions           // channel: [ versions.yml ]
}

