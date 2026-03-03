include { PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES } from '../../../modules/nf-core/pbcpgtools/alignedbamtocpgscores/main'
include { METHBAT_PROFILE                  } from '../../../modules/nf-core/methbat/profile/main'
include { CSVTK_MUTATE2                    } from '../../../modules/nf-core/csvtk/mutate2/main'
include { CSVTK_CONCAT                     } from '../../../modules/nf-core/csvtk/concat/main'
include { CSVTK_SORT                       } from '../../../modules/nf-core/csvtk/sort/main'

workflow CALL_METHYLATION_METHBAT {
    take:
    ch_bam_bai // channel: [ val(meta), bam, bai ]
    ch_regions // channel: [ val(meta), tsv ]

    main:
    ch_versions = channel.empty()

    PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES(
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.versions)

    PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.combined_bed
        .mix(
            PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.combined_bed_index,
            PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap1_bed,
            PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap1_bed_index,
            PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap2_bed,
            PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap2_bed_index,
        )
        .groupTuple()
        .set { ch_methbat_profile_in }

    METHBAT_PROFILE(
        ch_methbat_profile_in,
        ch_regions
    )
    ch_versions = ch_versions.mix(METHBAT_PROFILE.out.versions)

    /*
     * The region_profile files do not contain a sample identifier.
     * Therefore, we add it so we can concatenate files from multiple samples together in the next steps.
     */
    // TODO: Make this a separate workflow
    CSVTK_MUTATE2(
        METHBAT_PROFILE.out.region_profile,
        'tsv',
        'tsv'
    )

    CSVTK_MUTATE2.out.output
        .map { meta, region_profile_with_sample_id ->
            [ [ id: meta.family_id ], region_profile_with_sample_id ]
        }
        .groupTuple()
        .set { region_profiles_with_sample_id_per_family }

    // TODO: Add ext.prefix
    CSVTK_CONCAT(
        region_profiles_with_sample_id_per_family,
        'tsv',
        'tsv',
    )

    // TODO: Add ext.prefix
    CSVTK_SORT(
        CSVTK_CONCAT.out.csv,
        'tsv',
        'tsv',
    )

    emit:
    region_profile        = METHBAT_PROFILE.out.region_profile                   // channel: [ val(meta), path(tsv) ]
    asm_bed               = METHBAT_PROFILE.out.asm_bed                          // channel: [ val(meta), path(bed) ]
    pbcpg_biwgig_combined = PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.combined_bigwig // channel: [ val(meta), path(combined.bw) ]
    pbcpg_biwgig_hap1     = PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap1_bigwig     // channel: [ val(meta), path(hap1.bw) ]
    pbcpg_biwgig_hap2     = PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap2_bigwig     // channel: [ val(meta), path(hap2.bw) ]
    versions              = ch_versions                                          // channel: [ versions.yml ]
}
