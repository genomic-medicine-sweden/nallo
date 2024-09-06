include { PARAPHASE } from '../../modules/nf-core/paraphase/main'

workflow CALL_PARALOGS {

    take:
    bam_bai // channel: [ val(meta), bam, bai ]
    fasta   // channel: [ val(meta), fasta ]

    main:
    ch_versions = Channel.empty()

    PARAPHASE ( bam_bai, fasta, [[],[]] )
    ch_versions = ch_versions.mix(PARAPHASE.out.versions)

    emit:
    bam      = PARAPHASE.out.bam           // channel: [ val(meta), path(bam) ]
    bai      = PARAPHASE.out.bai           // channel: [ val(meta), path(bai) ]
    json     = PARAPHASE.out.json          // channel: [ val(meta), path(json) ]
    vcf      = PARAPHASE.out.vcf           // channel: [ val(meta), path(vcfs) ]
    tbi      = PARAPHASE.out.vcf_index     // channel: [ val(meta), path(tbis) ]
    versions = ch_versions                 // channel: [ versions.yml ]
}

