include { PARAPHASE        } from '../../modules/nf-core/paraphase/main'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow CALL_PARALOGS {

    take:
    bam_bai // channel: [ val(meta), bam ]
    fasta   // channel: [ val(meta), fasta ]

    main:
    ch_versions = Channel.empty()

    PARAPHASE ( bam_bai, fasta, [[],[]] )
    ch_versions = ch_versions.mix(PARAPHASE.out.versions)

    PARAPHASE.out.vcf
        .transpose() // Does create ~160 jobs per sample by default in hg38
        .set { bgzip_paraphase_vcfs }

    TABIX_BGZIPTABIX ( bgzip_paraphase_vcfs )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
    bam      = PARAPHASE.out.bam           // channel: [ val(meta), bam ]
    bai      = PARAPHASE.out.bai           // channel: [ val(meta), bai ]
    json     = PARAPHASE.out.json          // channel: [ val(meta), json ]
    vcf      = TABIX_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), gz, tbi ]

    versions = ch_versions                 // channel: [ versions.yml ]
}

