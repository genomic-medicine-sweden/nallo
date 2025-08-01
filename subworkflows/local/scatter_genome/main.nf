include { BEDTOOLS_MERGE           } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT            } from '../../../modules/nf-core/bedtools/sort/main'
include { BUILD_INTERVALS          } from '../../../modules/local/build_intervals/main'
include { BEDTOOLS_SPLIT           } from '../../../modules/nf-core/bedtools/split/main'

workflow SCATTER_GENOME {

    take:
    ch_fai             // channel: [optional] [ val(meta), path(fai) ]
    ch_input_bed       // channel: [optional] [ val(meta), path(bed) ]
    make_bed_from_fai  //    bool: Should we build a bed file from the fai?
    make_bed_intervals //    bool: Should we create intervals from the bed file?
    split_n            // integer: split bed into n regions

    main:
    ch_versions = Channel.empty()
    ch_bed = Channel.empty()
    ch_bed_intervals = Channel.empty()

    //
    // If no BED-file is provided then build intervals from reference
    //
    if( make_bed_from_fai ) {

        BUILD_INTERVALS (
            ch_fai
        )
        ch_versions = ch_versions.mix(BUILD_INTERVALS.out.versions)

        BUILD_INTERVALS.out.bed
            .set{ ch_bed }
    } else {
        ch_input_bed
            .set{ ch_bed }
    }

    //
    // Merge overlapping and then split BED regions for SNV calling
    //
    if( make_bed_intervals ) {

        if( split_n < 1 ) { error "Can't split bed file into less than one file" }

        // Sort and merge overlapping regions
        BEDTOOLS_SORT (
            ch_bed,
            []
        )
        ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

        BEDTOOLS_MERGE (
            BEDTOOLS_SORT.out.sorted
        )
        ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

        BEDTOOLS_SPLIT(
            BEDTOOLS_MERGE.out.bed.map { meta, bed ->
                [ meta, bed, split_n ]
            }
        )
        ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions)

        // Create a channel with the bed file and the total number of intervals (for groupKey)
        BEDTOOLS_SPLIT.out.beds
            .map { _meta, beds -> beds }
            .collect()
            .map{ it -> [ it, it.size() ] }
            .transpose()
            .set { ch_bed_intervals }

        // Since we don't check beforehand how many intervals it's possible to split the bed file into,
        // it could be that the number of intervals is less than the requested split_n.
        // This can happen if the bed file has too few regions.
        // We check this here, so it doens't fail later in the pipeline.
        ch_bed_intervals
            .count()
            .map { count ->
                if (count != split_n) {
                    error "Expected ${split_n}, but got ${count} files from splitting the BED file. " +
                          "This can happen if the input BED file (or fasta file if not using BED) has too few regions to split into ${split_n} parts. " +
                          "Please check the input files or set `--snv_calling_processes` to ${count}."
                }
            }
    }

    emit:
    bed           = ch_bed            // channel: [ val(meta), path(bed) ]
    bed_intervals = ch_bed_intervals  // channel: [ path(bed), val(num_intervals) ]
    versions      = ch_versions       // channel: [ versions.yml ]
}
