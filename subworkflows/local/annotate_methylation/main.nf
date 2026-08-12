include { CSVTK_CONCAT  } from '../../../modules/nf-core/csvtk/concat/main'
include { CSVTK_JOIN    } from '../../../modules/nf-core/csvtk/join/main'
include { CSVTK_MUTATE2 } from '../../../modules/nf-core/csvtk/mutate2/main'
include { CSVTK_SORT    } from '../../../modules/nf-core/csvtk/sort/main'

workflow ANNOTATE_METHYLATION {
    take:
    ch_region_profile // channel: [ val(meta), path(tsv) ]
    ch_map // channel: [ val(meta), path(tsv) ]

    main:
    CSVTK_JOIN(
        ch_region_profile.combine(ch_map).map { meta, region_profile, _meta_map, map ->
            [meta, [region_profile, map]]
        }
    )

    CSVTK_MUTATE2(
        CSVTK_JOIN.out.out_file,
        'tsv',
        'tsv',
    )

    ch_region_profiles_with_sample_id_per_family = CSVTK_MUTATE2.out.output
        .map { meta, region_profile_with_sample_id ->
            [[id: meta.family_id], region_profile_with_sample_id]
        }
        .groupTuple()

    CSVTK_CONCAT(
        ch_region_profiles_with_sample_id_per_family,
        'tsv',
        'tsv',
    )

    CSVTK_SORT(
        CSVTK_CONCAT.out.csv,
        'tsv',
        'tsv',
    )

    emit:
    methylation_annotation = CSVTK_SORT.out.sorted // channel: [ val(meta), path(tsv) ]
}
