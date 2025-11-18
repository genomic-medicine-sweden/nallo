include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'
include { GENERATE_GENS_DATA } from '../../../modules/local/generate_gens_data/main'
include { PREPROCESS_GENS_COV_INPUT } from '../../../modules/local/preprocess_gens_cov_input/main'

include { GATK4_DENOISEREADCOUNTS } from '../../../modules/nf-core/gatk4/denoisereadcounts/main'

// Thinking point: Do we want upd? meta? here

workflow PREPARE_GENS_INPUTS {
    take:
    ch_bam              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_gvcf             // channel: [mandatory] [ val(meta), path(gvcf), path(gvcf_tbi), path(baf_positions) ]
    gatk_header         // channel: [mandatory] [ path(txt) ]
    panel_of_normals    // channel: [mandatory] [ path(hd5) ]

    main:
    ch_versions = channel.empty()

    // meta, bam, bai, bed
    ch_mosdepth_in = ch_bam.map { meta, bam, bai -> tuple(meta, bam, bai, []) }
    ch_empty_fasta = ch_bam.map { meta, _bam, _bai -> tuple(meta, []) }

    // FIXME: Some questions to Viktor to be asked
    // Configure mosdepth to run with flags
    // -n (?)
    // --fast-mode (?)
    // -t 4 (?)
    // --by 100 (yes)
    // $SAMPLE_ID ?

    ch_mosdepth_in.view()
    ch_empty_fasta.view()

    // To configure
    // mosdepth -n --fast-mode -t 4 --by 100 test hg002_chr21_subsamp01_subset.bam
    MOSDEPTH(ch_mosdepth_in, ch_empty_fasta)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    MOSDEPTH.out.regions_bed.view()

    // What output?
    // Some proper configuration needed here

    channel
        .fromPath(panel_of_normals)
        .map { f -> tuple([], f) }
        .set { ch_pon }

    GATK4_DENOISEREADCOUNTS(
        MOSDEPTH.out.regions_bed,
        ch_pon
        // meta, pon
    )
    ch_versions = ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions)

    // PlotDenoisedCopyRatios
    // Cleanup

    // Consider using GATK4_DENOISEREADCOUNTS.out.standardized as well here
    // Also worth considering showing the raw data or normalized raw
    // without the PoN
    PREPROCESS_GENS_COV_INPUT(
        GATK4_DENOISEREADCOUNTS.out.standardized,
        gatk_header,
    )
    ch_versions = ch_versions.mix(PREPROCESS_GENS_COV_INPUT.out.versions)

    ch_gens_input = PREPROCESS_GENS_COV_INPUT.out.output.join(ch_gvcf)

    // Input: meta, coverage, gvcf
    GENERATE_GENS_DATA(
        ch_gens_input
    )
    ch_versions = ch_versions.mix(GENERATE_GENS_DATA.out.versions)

    // mosdepth
    // awk to GATK format
    // Python script

    // FIXME: Do we want a "without normal" cov output as well?
    emit:
    cov_bed_tbi = GENERATE_GENS_DATA.out.cov_bed_tbi
    baf_bed_tbi = GENERATE_GENS_DATA.out.baf_bed_tbi
    versions = ch_versions
}