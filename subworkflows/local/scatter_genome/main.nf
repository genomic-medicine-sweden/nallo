include { BEDTOOLS_MERGE               } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_SORT                } from '../../../modules/nf-core/bedtools/sort/main'
include { GAWK as GAWK_BUILD_INTERVALS } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_SPLIT               } from '../../../modules/nf-core/bedtools/split/main'
include { GAWK as GAWK_EXTRACT_REGIONS } from '../../../modules/nf-core/gawk/main'

workflow SCATTER_GENOME {
    take:
    ch_fai             // channel: [optional] [ val(meta), path(fai) ]
    ch_input_bed       // channel: [optional] [ val(meta), path(bed) ]
    make_bed_from_fai  //    bool: Should we build a bed file from the fai?
    split_n            // integer: split bed into n regions

    main:
    ch_versions = channel.empty()
    ch_bed = channel.empty()
    ch_bed_intervals = channel.empty()
    ch_bed_nuclear_mitochondrial_intervals = channel.empty()

    /*
    * If make_bed_from_fai is true then build intervals from reference
    */
    if (make_bed_from_fai) {

        GAWK_BUILD_INTERVALS(
            ch_fai,
            [],
            false,
        )

        GAWK_BUILD_INTERVALS.out.output.set { ch_bed }
    }
    else {
        ch_input_bed.set { ch_bed }
    }

    // Sort and merge overlapping regions
    BEDTOOLS_SORT(
        ch_bed,
        [],
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    BEDTOOLS_MERGE(
        BEDTOOLS_SORT.out.sorted
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

    // Add meta.genome so we can extract the mitochondrial region from the BED file
    BEDTOOLS_MERGE.out.bed
        .flatMap { meta, bed ->
            [[meta + [genome: "nuclear"], bed], [meta + [genome: "mitochondrial"], bed]]
        }
        .set { ch_input_gawk }

    // Exctract according to meta.genome, logic is in the config file
    GAWK_EXTRACT_REGIONS(
        ch_input_gawk,
        [],
        false,
    )

    GAWK_EXTRACT_REGIONS.out.output
        .branch { meta, _bed ->
            mitochondrial: meta.genome == "mitochondrial"
            nuclear: meta.genome == "nuclear"
        }
        .set { ch_bed_genomes }

    add_bed_count(ch_bed_genomes.nuclear)
        .map { meta, bed, num_intervals -> [meta.subMap('genome'), bed, num_intervals] }
        .set { ch_bed_intervals }

    add_bed_count(ch_bed_genomes.nuclear
        .mix(ch_bed_genomes.mitochondrial))
        .map { meta, bed, num_intervals -> [meta.subMap('genome'), bed, num_intervals] }
        .set { ch_bed_nuclear_mitochondrial_intervals }

    // Make bed interval if split_n > 1, otherwise just pass the bed file through
    if (split_n > 1) {

        // Split the nuclear bed file into n regions for SNV calling
        BEDTOOLS_SPLIT(
            ch_bed_genomes.nuclear.map { meta, bed ->
                [meta, bed, split_n]
            }
        )
        ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions)

        // Transpose the output so that we have [ val(meta), path(bed), val(num_intervals) ] for each interval file (chunk)
        BEDTOOLS_SPLIT.out.beds
            .transpose()
            .set { ch_bed_intervals }

        // Add the bed count and mix the nuclear intervals and mitochondrial channels
        add_bed_count(ch_bed_intervals
            .mix(ch_bed_genomes.mitochondrial))
            .map { meta, bed, num_intervals -> [meta.subMap('genome'), bed, num_intervals] }
            .set { ch_bed_nuclear_mitochondrial_intervals }

        /*
         * Since we don't check beforehand how many intervals it's possible to split the bed file into,
         * it could be that the number of intervals is less than the requested split_n.
         * This can happen if the bed file has too few regions.
         * We check this here, so it doesn't fail later in the pipeline.
         */
        ch_bed_intervals
            .count()
            .map { count ->
                if (count != split_n) {
                    error(
                        "Expected ${split_n}, but got ${count} files from splitting the BED file. " + "This can happen if the input BED file (or fasta file if not using BED) has too few regions to split into ${split_n} parts. " + "Please check the input files or set `--snv_calling_processes` to ${count}."
                    )
                }
            }
    }
    else if (split_n < 1) {
        error("Cannot split bed into ${split_n} regions, split_n should be minimum 1.")
    }

    emit:
    bed                                 = BEDTOOLS_MERGE.out.bed                 // channel: [ val(meta), path(bed) ]
    bed_nuclear_intervals               = ch_bed_intervals                       // channel: [ val(meta), path(bed), val(num_intervals) ]
    bed_nuclear_mitochondrial_intervals = ch_bed_nuclear_mitochondrial_intervals // channel: [ val(meta), path(bed), val(num_intervals) ]
    versions                            = ch_versions                            // channel: [ versions.yml ]
}
// Function to add the bed count to a channel
def add_bed_count(channel_with_beds) {
    def bed_count = channel_with_beds
        .map { _meta, bed -> bed }
        .collect()
        .map { beds -> beds.size() }

    channel_with_beds.combine(bed_count)
}
