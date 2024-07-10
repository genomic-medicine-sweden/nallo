include { BEDTOOLS_MERGE           } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT            } from '../../modules/nf-core/bedtools/sort/main'
include { BUILD_INTERVALS          } from '../../modules/local/build_intervals/main'
include { MINIMAP2_INDEX           } from '../../modules/nf-core/minimap2/index/main'
include { GUNZIP as GUNZIP_FASTA   } from '../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX           } from '../../modules/nf-core/samtools/faidx/main'
include { SPLIT_BED_CHUNKS         } from '../../modules/local/split_bed_chunks/main'
include { UNTAR as UNTAR_VEP_CACHE } from '../../modules/nf-core/untar/main'

workflow PREPARE_GENOME {

    take:
    fasta_in     // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_vep_cache // channel: [optional] [ path(cache) ]
    ch_input_bed // channel: [optional] [ val(meta, path(bed) ]

    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.empty()
    ch_bed = Channel.empty()

    fasta_file = fasta_in.map{meta, file -> file}

    // Will not catch cases where fasta is bgzipped
    if ( params.fasta.endsWith('.gz') ) {
        GUNZIP_FASTA(fasta_in)
            .gunzip
            .collect()
            .set{ch_fasta}

        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions.first())
    } else {
        fasta_in
            .set{ch_fasta}
    }

    SAMTOOLS_FAIDX ( ch_fasta, [[],[]] )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    MINIMAP2_INDEX ( ch_fasta )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    UNTAR_VEP_CACHE (ch_vep_cache)
    ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

    UNTAR_VEP_CACHE.out.untar
        .map { meta, files -> [files] }
        .collect()
        .set { untarred_vep }

    //
    // If no BED-file is provided then build intervals from reference
    //
    if( !params.bed ) {
        SAMTOOLS_FAIDX.out.fai
            .map{ name, fai -> [ ['id':name ], fai ] }
            .set{ ch_build_intervals_in }

        BUILD_INTERVALS( ch_build_intervals_in )
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
    if(!params.skip_short_variant_calling) {

        // Sort and merge overlapping regions
        BEDTOOLS_SORT ( ch_bed, [] )
        ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

        BEDTOOLS_MERGE ( BEDTOOLS_SORT.out.sorted )
        ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

        SPLIT_BED_CHUNKS( BEDTOOLS_MERGE.out.bed, params.parallel_snv )
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
    bed_intervals = ch_bed_intervals                   // channek: [ path(bed), val(num_intervals) ]
    fasta         = ch_fasta                           // channel: [ val(meta), path(fasta) ]
    fai           = SAMTOOLS_FAIDX.out.fai.collect()   // channel: [ val(meta), path(fai) ]
    mmi           = MINIMAP2_INDEX.out.index.collect() // channel: [ val(meta), path(mmi) ]
    vep_resources = untarred_vep                       // channel: [ path(cache) ]
    versions      = ch_versions                        // channel: [ versions.yml ]
}
