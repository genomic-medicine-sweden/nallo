include { ADD_FOUND_IN_TAG             } from '../../modules/local/add_found_in_tag'
include { HIFICNV                      } from '../../modules/local/pacbio/hificnv'

workflow CNV {

    take:
    ch_bam_bai_vcf // channel: [ val(meta), [[ bam ], [bai], [vcf] ]
    ch_fasta
    ch_expected_xy_bed
    ch_expected_xx_bed
    ch_exclude_bed

    main:
    ch_versions     = Channel.empty()

    // Split samples if male or female
    ch_bam_bai_vcf
        .branch{ meta, bam, bai, vcf ->
            male:   meta.sex == 1
            female: meta.sex == 2
        }
        .set{branched_bam}

    // Then merge male with XY-bed
    branched_bam
        .male
        .combine(ch_expected_xy_bed)
        .set {branched_bam_male}

    // And female with XX-bed
    branched_bam
        .female
        .combine(ch_expected_xx_bed)
        .set {branched_bam_female}

    // Concatenate the two channels
    branched_bam_male
        .concat(branched_bam_female)
        .set{ ch_hificnv_in }

    // Run HiFiCNV
    HIFICNV(ch_hificnv_in, ch_fasta, ch_exclude_bed)
    ch_versions = ch_versions.mix(HIFICNV.out.versions)

    // Add FOUND_IN=hificnv to VCF
    ADD_FOUND_IN_TAG (
        HIFICNV.out.vcf.map { meta, vcf -> [ meta, vcf, [] ] },
        "hificnv"
    )
    ch_versions = ch_versions.mix(ADD_FOUND_IN_TAG.out.versions)

    emit:
    vcf = ADD_FOUND_IN_TAG.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi = ADD_FOUND_IN_TAG.out.tbi // channel: [ val(meta), path(tbi) ]
    versions = ch_versions         // channel: [ versions.yml ]
}

