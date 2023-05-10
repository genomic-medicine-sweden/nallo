// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { HIFIASM } from '../../modules/local/hifiasm.nf'
include { GFA2FASTA } from '../../modules/local/gfa2fasta.nf'

workflow ASSEMBLY {

    take:
    // TODO nf-core: edit input (take) channels
    ch_reads // channel: [ val(meta), [ bam ] ]

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

