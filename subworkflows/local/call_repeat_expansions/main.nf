include { ADD_FOUND_IN_TAG                         } from '../../../modules/local/add_found_in_tag/main'
include { TRGT_GENOTYPE                            } from '../../../modules/nf-core/trgt/genotype/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRGT    } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRGT      } from '../../../modules/nf-core/samtools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_TRGT      } from '../../../modules/nf-core/bcftools/sort/main'
include { TRGT_MERGE                               } from '../../../modules/nf-core/trgt/merge/main'
include { BCFTOOLS_INDEX                           } from '../../../modules/nf-core/bcftools/index/main'
include { STRDUST                                  } from '../../../modules/nf-core/strdust/'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_STRDUST } from '../../../modules/bcftools/merge/'

workflow CALL_REPEAT_EXPANSIONS {

    take:
    ch_bam_bai  // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai      // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed // channel: [mandatory] [ val(meta), path(bed) ]
    str_caller  // string

    main:
    ch_versions         = Channel.empty()
    ch_sample_vcf       = Channel.empty()
    ch_sample_tbi       = Channel.empty()
    ch_family_vcf       = Channel.empty()
    ch_family_tbi       = Channel.empty()
    ch_sample_bam       = Channel.empty()
    ch_sample_bai       = Channel.empty()

    if (str_caller == "trgt") {
        ch_bam_bai
            .map { meta, bam, bai -> [ meta, bam, bai, meta.sex == 1 ? 'XY' : 'XX' ] }
            .set { ch_trgt_input }

        // Run TRGT
        TRGT_GENOTYPE (
            ch_trgt_input,
            ch_fasta,
            ch_fai,
            ch_bed
        )
        ch_versions = ch_versions.mix(TRGT_GENOTYPE.out.versions)

        // Sort and index bam
        SAMTOOLS_SORT_TRGT (
            TRGT_GENOTYPE.out.bam,
            [[],[]]
        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_TRGT.out.versions)

        SAMTOOLS_INDEX_TRGT ( SAMTOOLS_SORT_TRGT.out.bam )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TRGT.out.versions)

        // Add FOUND_IN=TRGT tag
        ADD_FOUND_IN_TAG (
            TRGT_GENOTYPE.out.vcf.map { meta, vcf -> [ meta, vcf, [] ] },
            "TRGT"
        )

        // Sort and index bcf
        BCFTOOLS_SORT_TRGT ( ADD_FOUND_IN_TAG.out.vcf )
        ch_versions = ch_versions.mix(BCFTOOLS_SORT_TRGT.out.versions)

        BCFTOOLS_SORT_TRGT.out.vcf
            .join( BCFTOOLS_SORT_TRGT.out.tbi, failOnMismatch:true, failOnDuplicate:true )
            .map { meta, bcf, csi -> [ [ id : meta.family_id ], bcf, csi ] }
            .groupTuple()
            .set{ ch_trgt_merge_in }

        TRGT_MERGE (
            ch_trgt_merge_in,
            [[],[]],
            [[],[]],
        )
        ch_versions = ch_versions.mix(TRGT_MERGE.out.versions)

        BCFTOOLS_INDEX(
            TRGT_MERGE.out.vcf
        )
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

        ch_sample_vcf  = ch_sample_vcf.mix(BCFTOOLS_SORT_TRGT.out.vcf)
        ch_sample_tbi  = ch_sample_tbi.mix(BCFTOOLS_SORT_TRGT.out.tbi)
        ch_family_vcf  = ch_family_vcf.mix(TRGT_MERGE.out.vcf)
        ch_family_tbi  = ch_family_tbi.mix(BCFTOOLS_INDEX.out.tbi)
        ch_sample_bam  = ch_sample_bam.mix(SAMTOOLS_SORT_TRGT.out.bam)
        ch_sample_bai  = ch_sample_bai.mix(SAMTOOLS_INDEX_TRGT.out.bai)

    } else if (str_caller == "strdust"){
        STRDUST (
            ch_bam_bai,
            ch_fasta,
            ch_fai,
            ch_bed
        )
        ch_versions.mix(STRDUST.out.versions)

        ADD_FOUND_IN_TAG (
            STRDUST.out.vcf.join(STRDUST.out.tbi),
            "STRdust"
        )

        ADD_FOUND_IN_TAG.out.vcf
            .join(ADD_FOUND_IN_TAG.out.tbi, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, vcf, tbi -> [ [ id: meta.family_id ], vcf, tbi ] }
            .groupTuple()
            .set { ch_bcftools_merge_in }

        BCFTOOLS_MERGE_STRDUST (
            ch_bcftools_merge_in,
            [ [], [] ],
            [ [], [] ],
            [ [], [] ]
        )
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)
        ch_sample_vcf  = ch_sample_vcf.mix(ADD_FOUND_IN_TAG.out.vcf)
        ch_sample_tbi  = ch_sample_tbi.mix(ADD_FOUND_IN_TAG.out.tbi)
        ch_family_vcf  = ch_family_vcf.mix(BCFTOOLS_MERGE.out.vcf)
        ch_family_tbi  = ch_family_tbi.mix(BCFTOOLS_INDEX.out.tbi)
    }

    emit:
    sample_vcf  = ch_sample_vcf // channel: [ val(meta), path(vcf) ]
    sample_tbi  = ch_sample_tbi // channel: [ val(meta), path(tbi) ]
    family_vcf  = ch_family_vcf // channel: [ val(meta), path(vcf) ]
    family_tbi  = ch_family_tbi // channel: [ val(meta), path(tbi) ]
    sample_bam  = ch_sample_bam // channel: [ val(meta), path(bam) ]
    sample_bai  = ch_sample_bai // channel: [ val(meta), path(bai) ]
    versions    = ch_versions   // channel: [ versions.yml ]
}

