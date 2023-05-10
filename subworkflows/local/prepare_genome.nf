// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {

    take:
    ch_fasta // channel: [ val(meta), fasta ]
    
    main:
    ch_versions = Channel.empty()

    SAMTOOLS_FAIDX ( ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    emit:
    fasta = ch_fasta                  // channel: [ val(meta), fasta ]
    fai   = SAMTOOLS_FAIDX.out.fai    // channel: [ val(meta), fai ]

    versions = ch_versions            // channel: [ versions.yml ]
}

