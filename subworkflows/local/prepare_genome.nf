// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {

    take:
    // TODO nf-core: edit input (take) channels
    ch_fasta // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SAMTOOLS_FAIDX ( ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    
    emit:
    // TODO nf-core: edit emitted channels
    fasta = ch_fasta
    fai      = SAMTOOLS_FAIDX.out.fai          // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

