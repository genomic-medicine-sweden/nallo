include { ADD_FOUND_IN_TAG } from '../../../modules/local/add_found_in_tag'
include { HIFICNV          } from '../../../modules/local/pacbio/hificnv'
include { SVDB_MERGE       } from '../../../modules/nf-core/svdb/merge/main'
include { TABIX_TABIX      } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_CNVS {

    take:
    ch_bam_bai_vcf     // channel: [ val(meta), path(bam), path(vcf) ]
    ch_fasta           // channel: [ val(meta), path(fasta) ]
    ch_expected_xy_bed // channel: [ val(meta), path(bed) ]
    ch_expected_xx_bed // channel: [ val(meta), path(bed) ]
    ch_exclude_bed     // channel: [ val(meta), path(bed) ]

    main:
    ch_versions     = Channel.empty()

    ch_bam_bai_vcf
        .map { meta, bam, bai, vcf -> [ meta, bam, bai, vcf, meta.sex ] }
        .set { ch_hificnv_in }

    // Run HiFiCNV
    HIFICNV (
        ch_hificnv_in,
        ch_fasta,
        ch_expected_xy_bed,
        ch_expected_xx_bed,
        ch_exclude_bed
    )
    ch_versions = ch_versions.mix(HIFICNV.out.versions)

    // Add FOUND_IN=hificnv to VCF
    ADD_FOUND_IN_TAG (
        HIFICNV.out.vcf.map { meta, vcf -> [ meta, vcf, [] ] },
        "hificnv"
    )
    ch_versions = ch_versions.mix(ADD_FOUND_IN_TAG.out.versions)

    ADD_FOUND_IN_TAG.out.vcf
        .map { meta, vcf -> [ [ 'id': meta.family_id ], vcf ] }
        .groupTuple()
        .set { svdb_merge_in }

    // Merge the files
    SVDB_MERGE (
        svdb_merge_in,
        [],
        true
    )
    ch_versions = ch_versions.mix(SVDB_MERGE.out.versions)

    TABIX_TABIX ( SVDB_MERGE.out.vcf )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    sample_vcf  = ADD_FOUND_IN_TAG.out.vcf // channel: [ val(meta), path(vcf) ]
    sample_tbi  = ADD_FOUND_IN_TAG.out.tbi // channel: [ val(meta), path(tbi) ]
    family_vcf  = SVDB_MERGE.out.vcf       // channel: [ val(meta), path(vcf) ]
    family_tbi  = TABIX_TABIX.out.tbi      // channel: [ val(meta), path(tbi) ]
    versions = ch_versions                 // channel: [ versions.yml ]
}

