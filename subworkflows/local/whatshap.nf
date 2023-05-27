include { WHATSHAP_PHASE } from '../../modules/local/whatshap/phase/main'
include { WHATSHAP_STATS } from '../../modules/local/whatshap/stats/main'
include { WHATSHAP_HAPLOTAG } from '../../modules/local/whatshap/haplotag/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'

workflow WHATSHAP {
    take:
        ch_vcf_bam_bai
        ch_fasta
        ch_fai
    main:
        ch_versions = Channel.empty()
        
        fasta = ch_fasta.map{it[1]}
        fai = ch_fai.map{it[1]}

        WHATSHAP_PHASE ( ch_vcf_bam_bai.combine(fasta).combine(fai) )
        ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions.first())
   
        WHATSHAP_STATS ( WHATSHAP_PHASE.out.vcf_tbi )
        ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions.first())
        
        WHATSHAP_PHASE.out.vcf_tbi
            .join(ch_vcf_bam_bai, by: 0)
            .map{ [ it[0], it[1], it[2], it[4], it[5] ] }
            .combine(fasta)
            .combine(fai)
            .set{ ch_whatshap_haplotag_in }
    
        WHATSHAP_HAPLOTAG( ch_whatshap_haplotag_in )
        ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions.first())


        SAMTOOLS_INDEX ( WHATSHAP_HAPLOTAG.out.bam )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        // Combine haplotagged bams with bai
        WHATSHAP_HAPLOTAG
            .out.bam
            .concat(SAMTOOLS_INDEX.out.bai)
            .groupTuple().flatten().collate(3)
            .set{ch_bam_bai_haplotagged}

    emit:
    haplotagged_bam_bai = ch_bam_bai_haplotagged
    versions = ch_versions
}
