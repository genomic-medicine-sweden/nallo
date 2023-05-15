include { DIPCALL } from '../../modules/local/dipcall'

workflow ASSEMBLY_VARIANT_CALLING {

    take:
    ch_haplotypes 
    ch_fasta
    ch_fai
    ch_mmi
    ch_ped
    ch_par

    main:
    ch_sv_calls_vcf = Channel.empty()
    ch_versions     = Channel.empty()
    
    dipcall_input = ch_haplotypes.flatten().collate(3).combine(ch_fasta.map{ it[1] }).combine( ch_fai.map{it[1]} ).combine(ch_mmi.map{it[1]}).combine(ch_par).combine(Channel.value(1))
    
    // Make sure reference has chrY PARs hard masked
    DIPCALL ( dipcall_input )
 
    ch_versions = ch_versions.mix(DIPCALL.out.versions)
    
    emit:
    ch_sv_calls_vcf        // channel: ? 
    versions = ch_versions // channel: [ versions.yml ]
}

