include { WHATSHAP_PHASE                            } from '../../modules/local/whatshap/phase/main'
include { WHATSHAP_STATS                            } from '../../modules/local/whatshap/stats/main'
include { WHATSHAP_HAPLOTAG                         } from '../../modules/local/whatshap/haplotag/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_WHATSHAP } from '../../modules/nf-core/samtools/index/main'
include { CRAMINO as CRAMINO_PHASED                 } from '../../modules/local/cramino'

workflow WHATSHAP {
    take:
        ch_vcf
        ch_bam_bai
        fasta
        fai
    main:
        ch_versions = Channel.empty()
    
        // Phase VCF
        WHATSHAP_PHASE ( ch_vcf.join(ch_bam_bai), fasta, fai )
        // Get phased stats 
        WHATSHAP_STATS ( WHATSHAP_PHASE.out.vcf_tbi )
         
        WHATSHAP_PHASE.out.vcf_tbi
            .join(ch_bam_bai)
            .set{ ch_whatshap_haplotag_in }
        
        // Haplotag reads
        WHATSHAP_HAPLOTAG(ch_whatshap_haplotag_in, fasta, fai)

        // Index reads
        SAMTOOLS_INDEX_WHATSHAP ( WHATSHAP_HAPLOTAG.out.bam )
        
        // Combine haplotagged bams with bai
        WHATSHAP_HAPLOTAG
            .out.bam
            .join(SAMTOOLS_INDEX_WHATSHAP.out.bai)
            .set{ch_bam_bai_haplotagged}
         
         // Prepare inputs
        ch_bam_bai_haplotagged
            .map{ meta, bam, bai -> [ meta, bam, bai, [] ] }
            .set{ ch_mosdepth_in }   
    
        CRAMINO_PHASED(ch_bam_bai_haplotagged)

        // Get versions
        ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions.first())
        ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions.first())
        ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_WHATSHAP.out.versions.first())
        ch_versions = ch_versions.mix(CRAMINO_PHASED.out.versions.first())

    emit:
    haplotagged_bam_bai = ch_bam_bai_haplotagged // channel: [ val(meta), bam, bai ]
    versions            = ch_versions            // channel: [ versions.yml ]
}
