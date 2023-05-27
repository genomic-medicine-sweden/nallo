include { SAMTOOLS_INDEX                                   } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX    as SAMTOOLS_INDEX_HP1          } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX    as SAMTOOLS_INDEX_HP2          } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW     as SAMTOOLS_VIEW_HP1           } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW     as SAMTOOLS_VIEW_HP2           } from '../../modules/nf-core/samtools/view/main'
include { MODKIT_PILEUP                                    } from '../../modules/local/modkit/pileup/main'
include { MODKIT_PILEUP     as MODKIT_PILEUP_HP1           } from '../../modules/local/modkit/pileup/main'
include { MODKIT_PILEUP     as MODKIT_PILEUP_HP2           } from '../../modules/local/modkit/pileup/main'

include { WHATSHAP } from '../../subworkflows/local/whatshap'

workflow METHYLATION {

    take:
    ch_bam_bai // channel: [ val(meta), [[ bam ], [bai]] ]
    ch_vcf
    ch_fasta
    ch_fai

    main:
    ch_versions           = Channel.empty()
    
    // TODO: Allow PED for more accurate phasing 
    // TODO: Move haplotagging to its own workflow?
    // Phase SNPs-calls - move to SV-calling subworkflow? 
    
    ch_vcf
        .join(ch_bam_bai, by: 0)
        .map { meta, vcf, bam, bai -> 
            return [meta, vcf, bam, bai]}
        .set{ ch_whatshap_phase_in }
    
    WHATSHAP ( ch_whatshap_phase_in, ch_fasta, ch_fai)

    fa = WHATSHAP.out.haplotagged_bam_bai.combine(ch_fasta.map{it[1]}).map{ [ it[0], it[3] ]}

    // Extract haplotagged reads into separate files 
    SAMTOOLS_VIEW_HP1( WHATSHAP.out.haplotagged_bam_bai, fa, [] )
    SAMTOOLS_VIEW_HP2( WHATSHAP.out.haplotagged_bam_bai, fa, [] )

    // Index new bams 
    SAMTOOLS_INDEX_HP1 ( SAMTOOLS_VIEW_HP1.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_HP1.out.versions.first())

    SAMTOOLS_INDEX_HP2 ( SAMTOOLS_VIEW_HP2.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_HP2.out.versions.first())
    
    // Combine HP1 bam bai 
    SAMTOOLS_VIEW_HP1
        .out.bam
        .concat(SAMTOOLS_INDEX_HP1.out.bai)
        .groupTuple().flatten().collate(3)
        .set{ch_bam_bai_hp1}
    
    // Combine HP2 bam bai
    SAMTOOLS_VIEW_HP2
        .out.bam
        .concat(SAMTOOLS_INDEX_HP2.out.bai)
        .groupTuple().flatten().collate(3)
        .set{ch_bam_bai_hp2}
     
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_HP1.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_HP2.out.versions.first())
    
    // Sometimes you may not want per haplotype methylation (e.g. chrX males)
    // ...Maybe we could to this only on male chrX and chrY?
    // TODO: Fix inconsistent naming between hp1/2 and haplotagged bams?
    MODKIT_PILEUP( WHATSHAP.out.haplotagged_bam_bai.combine(ch_fasta.map{it[1]}.combine(ch_fai.map{it[1]})))
    // Sometimes you may want per haplotype methylation
    MODKIT_PILEUP_HP1 ( ch_bam_bai_hp1.combine(ch_fasta.map{it[1]}.combine(ch_fai.map{it[1]})))
    MODKIT_PILEUP_HP2 ( ch_bam_bai_hp2.combine(ch_fasta.map{it[1]}.combine(ch_fai.map{it[1]})))
    
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)
    
    emit:
    ch_combined_bedmethyl = MODKIT_PILEUP.out.bed
    ch_combined_bedmethyl = MODKIT_PILEUP.out.bed
    emit:
    versions              = ch_versions                  // channel: [ versions.yml ]
}

