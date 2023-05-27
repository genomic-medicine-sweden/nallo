include { WHATSHAP_PHASE                                   } from '../../modules/local/whatshap/phase/main'
include { WHATSHAP_STATS                                   } from '../../modules/local/whatshap/stats/main'
include { WHATSHAP_HAPLOTAG                                } from '../../modules/local/whatshap/haplotag/main'
include { SAMTOOLS_INDEX                                   } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX    as SAMTOOLS_INDEX_HP1          } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX    as SAMTOOLS_INDEX_HP2          } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW     as SAMTOOLS_VIEW_HP1           } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW     as SAMTOOLS_VIEW_HP2           } from '../../modules/nf-core/samtools/view/main'
include { MODKIT_PILEUP                                    } from '../../modules/local/modkit/pileup/main'
include { MODKIT_PILEUP     as MODKIT_PILEUP_PER_HAPLOTYPE } from '../../modules/local/modkit/pileup/main'

workflow METHYLATION {

    take:
    ch_bam_bai // channel: [ val(meta), [[ bam ], [bai]] ]
    ch_vcf
    ch_fasta
    ch_fai

    main:
    ch_versions           = Channel.empty()
    ch_phased_vcf         = Channel.empty()
    ch_haplotagged_reads  = Channel.empty()
    ch_combined_bedmethyl = Channel.empty()
    ch_per_hap_bedmethyl  = Channel.empty()
    
    meta           = ch_bam_bai.map{ it[0] }
    fasta_fai      = ch_fasta.combine(meta).map{[it[2], it[1]]}.combine(ch_fai.map{ it[1] })
    
    // TODO: Allow PED for more accurate phasing 
    // TODO: Move haplotagging to its own workflow?
    // Phase SNPs-calls - move to SV-calling subworkflow? 
    WHATSHAP_PHASE ( ch_vcf, ch_bam_bai, fasta_fai)
    ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions.first())
    
    // Get stats
    WHATSHAP_STATS ( WHATSHAP_PHASE.out.vcf_tbi )
    ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions.first())
    
    // Haplotag reads
    WHATSHAP_HAPLOTAG( WHATSHAP_PHASE.out.vcf_tbi, ch_bam_bai, fasta_fai )
    ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions.first())
    
    // Index haplotagged bams 
    SAMTOOLS_INDEX ( WHATSHAP_HAPLOTAG.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Combine haplotagged bams with bai
    WHATSHAP_HAPLOTAG
        .out.bam
        .concat(SAMTOOLS_INDEX.out.bai)
        .groupTuple().flatten().collate(3)
        .set{ch_bam_bai_haplotagged}

    // Extract haplotagged reads into separate files 
    SAMTOOLS_VIEW_HP1( ch_bam_bai_haplotagged, fasta_fai.map{ [it[0], it[1]] }, [] )
    SAMTOOLS_VIEW_HP2( ch_bam_bai_haplotagged, fasta_fai.map{ [it[0], it[1]] }, [] )

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
    
    // Make double meta's for single HP-bams
    modkit_meta = ch_bam_bai_hp1.concat(ch_bam_bai_hp2).map{ it[0] }
    modkit_fasta_fai = modkit_meta.combine(fasta_fai.map{ [it[1], it[2]]})

    // Sometimes you may not want per haplotype methylation (e.g. chrX males)
    // ...Maybe we could to this only on male chrX and chrY?
    // TODO: Fix inconsistent naming between hp1/2 and haplotagged bams?
    MODKIT_PILEUP( ch_bam_bai_haplotagged, fasta_fai)

    // Sometimes you may want per haplotype methylation
    MODKIT_PILEUP_PER_HAPLOTYPE ( ch_bam_bai_hp1.concat(ch_bam_bai_hp2) , modkit_fasta_fai)
    
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)
    ch_versions = ch_versions.mix(MODKIT_PILEUP_PER_HAPLOTYPE.out.versions)
    
    emit:
    ch_phased_vcf         = WHATSHAP_PHASE.out.vcf_tbi
    ch_haplotagged_reads  = ch_bam_bai_haplotagged
    ch_combined_bedmethyl = MODKIT_PILEUP.out.bed
    ch_combined_bedmethyl = MODKIT_PILEUP.out.bed
    versions              = ch_versions                  // channel: [ versions.yml ]
}

