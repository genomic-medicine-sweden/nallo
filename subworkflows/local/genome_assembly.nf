include { HIFIASM   } from '../../modules/nf-core/hifiasm/main'
include { GFA2FASTA } from '../../modules/local/gfa2fasta.nf'

workflow ASSEMBLY {

    take:
    ch_reads // channel: [ val(meta), fastq ]

    main:

    ch_versions = Channel.empty()
    ch_phased_contigs_fa = Channel.empty()
    

    HIFIASM ( ch_reads ) // [ [meta], fastq ]

    // Split haplotypes, keep id   
    gfa_haplotypes = HIFIASM // [ [meta], hap_* ]
        .out
        .phased_contigs
        .map{ [ it[0], it[1], it[0], it[2] ] }
        .flatten()
        .collate(2)

    gfa_haplotypes.view()


    GFA2FASTA ( gfa_haplotypes ) 

    // Collect haplotype pairs per sample
    ch_dual_assembly_fa = GFA2FASTA // [ [meta], hap_1, hap_2 ] ]
        .out
        .phased_contigs_fa
        .groupTuple()

    
    ch_versions = ch_versions.mix(HIFIASM.out.versions.first())
    ch_versions = ch_versions.mix(GFA2FASTA.out.versions.first())
    
    
    emit:
    assembled_haplotypes = ch_dual_assembly_fa

    versions = ch_versions                     // channel: [ versions.yml ]
}

