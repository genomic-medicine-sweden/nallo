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
    
    // This should work both with and without --trio, but really should skip PED-file and just rely on samplesheet
    dipcall_input = ch_haplotypes
                    .flatten()
                    .collate(3)
                    .map{ [it[0]['id'], it[0], it[1], it[2] ]}
                    .join(ch_ped.map{ [ it['id'], it['sex'] ] }, by: 0)
                    .map{ [it[1], it[2], it[3], it[4]]}
    
    par = ch_par.collect()

    //Make sure reference has chrY PARs hard masked
    DIPCALL ( dipcall_input, ch_fasta, ch_fai, ch_mmi, par )
    ch_versions = ch_versions.mix(DIPCALL.out.versions)
    
    emit:
    ch_sv_calls_vcf        // channel: ? 
    versions = ch_versions // channel: [ versions.yml ]
}

