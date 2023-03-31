// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { PBMM2_ALIGN } from '../../modules/local/pbmm2/align'

workflow ALIGN_READS {

    take:
    // TODO nf-core: edit input (take) channels
    ch_fastq // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    PBMM2_ALIGN ( ch_fastq.combine(ch_fasta.map { it[1] }) )
    ch_versions = ch_versions.mix(PBMM2_ALIGN.out.versions.first())

    
    emit:
    // TODO nf-core: edit emitted channels
    bam_bai = PBMM2_ALIGN.out.bam_bai

    versions = ch_versions                     // channel: [ versions.yml ]
}

