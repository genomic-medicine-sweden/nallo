include { BCFTOOLS_MERGE   } from '../../../modules/nf-core/bcftools/merge/'
include { STRDUST          } from '../../../modules/nf-core/strdust/'
include { TABIX_TABIX      } from '../../../modules/nf-core/tabix/tabix/main'
include { VCFEXPRESS       } from '../../../modules/nf-core/vcfexpress/main'

workflow CALL_REPEAT_EXPANSIONS_STRDUST {

    take:
    ch_bam_bai              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta                // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai                  // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed                  // channel: [mandatory] [ val(meta), path(bed) ]
    ch_vcfexpress_prelude   // channel: [mandatory] [ val(meta), path(lua) ]

    main:
    ch_versions = channel.empty()

    STRDUST (
        ch_bam_bai,
        ch_fasta,
        ch_fai,
        ch_bed
    )
    ch_versions.mix(STRDUST.out.versions)

    ch_lua_file = ch_vcfexpress_prelude.map { meta, lua -> lua }

    VCFEXPRESS (
        STRDUST.out.vcf,
        ch_lua_file
    )

    TABIX_TABIX (
        VCFEXPRESS.out.vcf
    )

    VCFEXPRESS.out.vcf
        .join(TABIX_TABIX.out.index, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi -> [ [ id: meta.family_id ], vcf, tbi ] }
        .groupTuple()
        .set { ch_bcftools_merge_in }

    BCFTOOLS_MERGE (
        ch_bcftools_merge_in,
        [ [], [] ],
        [ [], [] ],
        [ [], [] ]
    )

    emit:
    sample_vcf  = STRDUST.out.vcf          // channel: [ val(meta), path(vcf) ]
    sample_tbi  = STRDUST.out.tbi          // channel: [ val(meta), path(tbi) ]
    family_vcf  = BCFTOOLS_MERGE.out.vcf   // channel: [ val(meta), path(vcf) ]
    family_tbi  = BCFTOOLS_MERGE.out.index // channel: [ val(meta), path(tbi) ]
    versions    = ch_versions              // channel: [ versions.yml ]

}
