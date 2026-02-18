include { BEDTOOLS_MERGE                       } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT                        } from '../../../modules/nf-core/bedtools/sort/main'
include { BUILD_INTERVALS                      } from '../../../modules/local/build_intervals/main'
include { BEDTOOLS_SPLIT                       } from '../../../modules/nf-core/bedtools/split/main'
include { GAWK as GAWK_EXTRACT_MT_REGIONS      } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_EXTRACT_NUCLEAR_REGIONS } from '../../../modules/nf-core/gawk/main'
workflow SCATTER_GENOME {

    take:
    ch_fai             // channel: [optional] [ val(meta), path(fai) ]
    ch_input_bed       // channel: [optional] [ val(meta), path(bed) ]
    make_bed_from_fai  //    bool: Should we build a bed file from the fai?
    make_bed_intervals //    bool: Should we create intervals from the bed file?
    split_n            // integer: split bed into n regions

    main:
    ch_versions = channel.empty()
    ch_bed = channel.empty()
    ch_bed_intervals = channel.empty()

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

        // Extract the mitochondrial region from BED before spliting into 40 regions
        GAWK_EXTRACT_MT_REGIONS (
            BEDTOOLS_MERGE.out.bed,
            [],
            false,
        )
        ch_versions = ch_versions.mix(GAWK_EXTRACT_MT_REGIONS.out.versions)
        GAWK_EXTRACT_MT_REGIONS.out.output.view()

        // Extract the nuclear genome regions from BED before spliting into 40 regions
        GAWK_EXTRACT_NUCLEAR_REGIONS (
            BEDTOOLS_MERGE.out.bed,
            [],
            false,
        )
        ch_versions = ch_versions.mix(GAWK_EXTRACT_NUCLEAR_REGIONS.out.versions)
        GAWK_EXTRACT_NUCLEAR_REGIONS.out.output.view()

        BEDTOOLS_SPLIT(
            GAWK_EXTRACT_NUCLEAR_REGIONS.out.output.map { meta, bed ->
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
        // We check this here, so it doesn't fail later in the pipeline.
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
    bed           = ch_bed                              // channel: [ val(meta), path(bed) ]
    bed_intervals = ch_bed_intervals                    // channel: [ path(bed), val(num_intervals) ]
    bed_mt        = GAWK_EXTRACT_MT_REGIONS.out.output  // channel: [ val(meta), path(bed) ]
    versions      = ch_versions                         // channel: [ versions.yml ]
}
