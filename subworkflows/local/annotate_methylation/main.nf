include { CSVTK_CONCAT                           } from '../../../modules/nf-core/csvtk/concat/main'
include { CSVTK_JOIN as CSVTK_JOIN_REGION_GENES  } from '../../../modules/nf-core/csvtk/join/main'
include { CSVTK_JOIN as CSVTK_JOIN_REGION_LABELS } from '../../../modules/nf-core/csvtk/join/main'
include { CSVTK_MUTATE2                          } from '../../../modules/nf-core/csvtk/mutate2/main'
include { CSVTK_SORT                             } from '../../../modules/nf-core/csvtk/sort/main'
include { GAWK as GAWK_REGION_GENES              } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_REGION_LABELS             } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_REGION_VALUES             } from '../../../modules/nf-core/gawk/main'

workflow ANNOTATE_METHYLATION {
    take:
    ch_methbat_profiles // channel: [optional]  [ val(meta), path(tsv) ]
    ch_modkit_dmr       // channel: [optional]  [ val(meta), path(bed) ]
    ch_regions          // channel: [mandatory] [ val(meta), path(tsv) ]
    ch_map              // channel: [mandatory] [ val(meta), path(tsv) ]

    main:
    //
    // Once per run: region labels, and the genes of each label from the map
    //
    GAWK_REGION_LABELS(
        ch_regions,
        [],
        false,
    )

    // The region file must come first: the program reads it as ARGV[1]
    GAWK_REGION_GENES(
        ch_regions.combine(ch_map).map { meta, regions, _meta_map, map -> [meta, [regions, map]] },
        [],
        false,
    )

    //
    // Per sample: values per region, spread to labels and genes
    //
    ch_region_values = ch_methbat_profiles
        .map { meta, profile -> [meta + [methylation_caller: 'methbat'], profile] }
        .mix(
            ch_modkit_dmr.map { meta, dmr -> [meta + [methylation_caller: 'modkit'], dmr] }
        )

    // Drops the region name column, which the tools set to chrom:start-end
    GAWK_REGION_VALUES(
        ch_region_values,
        [],
        false,
    )

    CSVTK_JOIN_REGION_LABELS(
        GAWK_REGION_VALUES.out.output.combine(GAWK_REGION_LABELS.out.output).map { meta, values, _meta_labels, labels -> [meta, [labels, values]] }
    )

    CSVTK_MUTATE2(
        CSVTK_JOIN_REGION_LABELS.out.out_file,
        'tsv',
        'tsv',
    )

    CSVTK_JOIN_REGION_GENES(
        CSVTK_MUTATE2.out.output.combine(GAWK_REGION_GENES.out.output).map { meta, values, _meta_genes, genes -> [meta, [values, genes]] }
    )

    ch_annotated_per_family = CSVTK_JOIN_REGION_GENES.out.out_file
        .map { meta, annotated ->
            [[id: meta.family_id, methylation_caller: meta.methylation_caller], annotated]
        }
        .groupTuple()

    CSVTK_CONCAT(
        ch_annotated_per_family,
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
