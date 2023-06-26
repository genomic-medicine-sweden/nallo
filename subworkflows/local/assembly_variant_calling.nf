include { DIPCALL                                  } from '../../modules/local/dipcall'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DIPCALL } from '../../modules/nf-core/samtools/index/main'
include { MINIMAP2_INDEX as MINIMAP2_INDEX_DIPCALL } from '../../modules/nf-core/minimap2/index/main'

workflow ASSEMBLY_VARIANT_CALLING {

    take:
    ch_haplotypes 
    ch_fasta
    ch_fai
    ch_par

    main:
    ch_sv_calls_vcf = Channel.empty()
    ch_versions     = Channel.empty()
    
    ch_haplotypes
        .map{ meta, hap1, hap2 -> [meta, hap1, hap2, meta.sex] }
        .set{ ch_dipcall_input }
    
    //Make sure reference has chrY PARs hard masked
    MINIMAP2_INDEX_DIPCALL(ch_fasta)
    
    DIPCALL( ch_dipcall_input, ch_fasta, ch_fai, MINIMAP2_INDEX_DIPCALL.out.index, ch_par )
    
    SAMTOOLS_INDEX_DIPCALL ( DIPCALL.out.bam.transpose() )
    
    ch_versions = ch_versions.mix(MINIMAP2_INDEX_DIPCALL.out.versions)
    ch_versions = ch_versions.mix(DIPCALL.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DIPCALL.out.versions)
    
    emit:
    ch_sv_calls_vcf        // channel: ? 
    versions = ch_versions // channel: [ versions.yml ]
}

