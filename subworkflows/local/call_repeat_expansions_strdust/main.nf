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
    ch_vcfexpress_prelude   // path: [mandatory] lua file

    main:

    STRDUST (
        ch_bam_bai,
        ch_fasta,
        ch_fai,
        ch_bed
    )

    VCFEXPRESS (
        STRDUST.out.vcf,
        ch_vcfexpress_prelude
    )

    TABIX_TABIX (
        VCFEXPRESS.out.vcf
    )

    VCFEXPRESS.out.vcf
        .join(TABIX_TABIX.out.index, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi -> [ [ id: meta.family_id ], vcf, tbi ] }
        .groupTuple()
        .map { meta, vcfs, tbis -> [ meta, vcfs, tbis, [] ] }
        .set { ch_bcftools_merge_in }

    BCFTOOLS_MERGE (
        ch_bcftools_merge_in,
        ch_fasta.join(ch_fai, failOnMismatch: true, failOnDuplicate: true).collect()
    )

    emit:
    sample_vcf  = VCFEXPRESS.out.vcf       // channel: [ val(meta), path(vcf) ]
    sample_tbi  = TABIX_TABIX.out.index    // channel: [ val(meta), path(tbi) ]
    family_vcf  = BCFTOOLS_MERGE.out.vcf   // channel: [ val(meta), path(vcf) ]
    family_tbi  = BCFTOOLS_MERGE.out.index // channel: [ val(meta), path(tbi) ]

}
