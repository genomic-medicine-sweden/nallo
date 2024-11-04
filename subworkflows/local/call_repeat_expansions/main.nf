include { ADD_FOUND_IN_TAG                       } from '../../../modules/local/add_found_in_tag/main'
include { TRGT                                   } from '../../../modules/local/trgt'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRGT  } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRGT    } from '../../../modules/nf-core/samtools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_TRGT    } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_MERGE                         } from '../../../modules/nf-core/bcftools/merge/main'

workflow CALL_REPEAT_EXPANSIONS {

    take:
    ch_bam_bai  // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai      // channel: [mandatory] [ val(meta), path(fai) ]
    ch_trgt_bed // channel: [mandatory] [ val(meta), path(bed) ]

    main:
    ch_repeat_calls_vcf = Channel.empty()
    ch_versions         = Channel.empty()

    ch_bam_bai
        .map { meta, bam, bai -> [meta, bam, bai, meta.sex] }
        .set { ch_trgt_input }

    // Run TGRT
    TRGT (
        ch_trgt_input,
        ch_fasta,
        ch_fai,
        ch_trgt_bed
    )
    ch_versions = ch_versions.mix(TRGT.out.versions)

    // Sort and index bam
    SAMTOOLS_SORT_TRGT (
        TRGT.out.bam,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_TRGT.out.versions)

    SAMTOOLS_INDEX_TRGT ( SAMTOOLS_SORT_TRGT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TRGT.out.versions)

    // Add FOUND_IN=TRGT tag
    ADD_FOUND_IN_TAG (
        TRGT.out.vcf.map { meta, vcf -> [ meta, vcf, [] ] },
        "TRGT"
    )

    // Sort and index bcf
    BCFTOOLS_SORT_TRGT ( ADD_FOUND_IN_TAG.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_TRGT.out.versions)

    BCFTOOLS_SORT_TRGT.out.vcf
        .join( BCFTOOLS_SORT_TRGT.out.tbi )
        .map { meta, bcf, csi -> [ [ id : meta.family_id ], bcf, csi ] }
        .groupTuple()
        .set{ ch_bcftools_merge_in }

    BCFTOOLS_MERGE (
        ch_bcftools_merge_in,
        ch_fasta,
        ch_fai,
        [[],[]]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    emit:
    sample_vcf  = BCFTOOLS_SORT_TRGT.out.vcf  // channel: [ val(meta), path(vcf) ]
    family_vcf  = BCFTOOLS_MERGE.out.vcf      // channel: [ val(meta), path(vcf) ]
    sample_bam  = SAMTOOLS_SORT_TRGT.out.bam  // channel: [ val(meta), path(bam) ]
    sample_bai  = SAMTOOLS_INDEX_TRGT.out.bai // channel: [ val(meta), path(bai) ]
    versions    = ch_versions                 // channel: [ versions.yml ]
}

