include { SAMTOOLS_INDEX                       } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_HP1 } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_HP2 } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW  as SAMTOOLS_VIEW_HP1  } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW  as SAMTOOLS_VIEW_HP2  } from '../../modules/nf-core/samtools/view/main'
include { MODKIT_PILEUP                        } from '../../modules/local/modkit/pileup/main'
include { MODKIT_PILEUP  as MODKIT_PILEUP_HP1  } from '../../modules/local/modkit/pileup/main'
include { MODKIT_PILEUP  as MODKIT_PILEUP_HP2  } from '../../modules/local/modkit/pileup/main'

workflow METHYLATION {

    take:
    ch_bam_bai // channel: [ val(meta), bam, bai ]
    ch_vcf     // channel: [ val(meta), vcf.gz ]
    ch_fasta   // channel: [ val(meta), fasta ]
    ch_fai     // channel: [ val(meta), fai ]

    main:
    ch_versions = Channel.empty()
    
    // TODO: Allow PED for more accurate phasing 
    // TODO: Move haplotagging to its own workflow?
    // Phase SNPs-calls - move to SV-calling subworkflow? 

    // Extract haplotagged reads into separate files 
    SAMTOOLS_VIEW_HP1( ch_bam_bai, ch_fasta, [] )
    SAMTOOLS_VIEW_HP2( ch_bam_bai, ch_fasta, [] )

    // Index new bams 
    SAMTOOLS_INDEX_HP1 ( SAMTOOLS_VIEW_HP1.out.bam )
    SAMTOOLS_INDEX_HP2 ( SAMTOOLS_VIEW_HP2.out.bam )
    
    // Combine HP1 bam bai 
    SAMTOOLS_VIEW_HP1
        .out.bam
        .join(SAMTOOLS_INDEX_HP1.out.bai)
        .set{ch_bam_bai_hp1}
    
    // Combine HP2 bam bai
    SAMTOOLS_VIEW_HP2
        .out.bam
        .join(SAMTOOLS_INDEX_HP2.out.bai)
        .set{ch_bam_bai_hp2}
    
    // Sometimes you may not want per haplotype methylation (e.g. chrX males)
    // ...Maybe we could to this only on male chrX and chrY?
    // TODO: Fix inconsistent naming between hp1/2 and haplotagged bams?
    
    MODKIT_PILEUP( ch_bam_bai, ch_fasta, ch_fai)
    
    // Sometimes you may want per haplotype methylation
    MODKIT_PILEUP_HP1 ( ch_bam_bai_hp1, ch_fasta, ch_fai)
    MODKIT_PILEUP_HP2 ( ch_bam_bai_hp2, ch_fasta, ch_fai)
    
    // Get versions
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_HP1.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_HP2.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_HP1.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_HP2.out.versions.first())
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)
    
    emit:
    ch_combined_bedmethyl = MODKIT_PILEUP.out.bed // channel: [ val(meta), bed.gz, tbi ]
    versions              = ch_versions           // channel: [ versions.yml ]
}

