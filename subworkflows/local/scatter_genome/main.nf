include { BEDTOOLS_MAKEWINDOWS     } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MERGE           } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT            } from '../../../modules/nf-core/bedtools/sort/main'
include { BUILD_INTERVALS          } from '../../../modules/local/build_intervals/main'
include { SPLIT_BED_CHUNKS         } from '../../../modules/local/split_bed_chunks/main'

workflow SCATTER_GENOME {

    take:
    ch_fai             // channel: [optional] [ val(meta), path(fai) ]
    ch_input_bed       // channel: [optional] [ val(meta), path(bed) ]
    make_bed_from_fai  // bool
    make_bed_intervals // bool
    split_n            // integer: split bed into n regions

    main:
    ch_versions = Channel.empty()
    ch_bed = Channel.empty()
    ch_bed_intervals = Channel.empty()

    //
    // If no BED-file is provided then build intervals from reference
    //
    if( make_bed_from_fai ) {


        BUILD_INTERVALS ( ch_fai.map { id, fai -> [ [ 'id': id ], fai ] } )
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
        BEDTOOLS_SORT ( ch_bed, [] )
        ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

        BEDTOOLS_MERGE ( BEDTOOLS_SORT.out.sorted )
        ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

        BEDTOOLS_MAKEWINDOWS ( BEDTOOLS_MERGE.out.bed )
        ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

        SPLIT_BED_CHUNKS( BEDTOOLS_MAKEWINDOWS.out.bed, split_n )
        ch_versions = ch_versions.mix(SPLIT_BED_CHUNKS.out.versions)

        // Create a channel with the bed file and the total number of intervals (for groupKey)
        SPLIT_BED_CHUNKS.out.split_beds
            .collect()
            .map{ it -> [ it, it.size() ] }
            .transpose()
            .set { ch_bed_intervals }
    }

    emit:
    bed           = ch_bed                             // channel: [ val(meta), path(bed) ]
    bed_intervals = ch_bed_intervals                   // channel: [ path(bed), val(num_intervals) ]
    versions      = ch_versions                        // channel: [ versions.yml ]
}

