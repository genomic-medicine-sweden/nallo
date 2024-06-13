include { PARAPHASE        } from '../../modules/nf-core/paraphase/main'
include { SAMTOOLS_INDEX   } from '../../modules/nf-core/samtools/index/main'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow CALL_PARALOGS {

    take:
    ch_bam   // channel: [ val(meta), bam ]
    ch_fasta // channel: [ val(meta), fasta ]

    main:

    ch_versions = Channel.empty()

    // Needs bai, not csi
    SAMTOOLS_INDEX ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam
        .join( SAMTOOLS_INDEX.out.bai )
        .set { paraphase_in }

    PARAPHASE ( paraphase_in, ch_fasta, [[],[]] )
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

