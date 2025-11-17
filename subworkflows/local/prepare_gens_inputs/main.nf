include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'
include { GENERATE_GENS_DATA } from '../../../modules/local/generate_gens_data/main'
include { PREPROCESS_GENS_COV_INPUT } from '../../../modules/local/preprocess_gens_cov_input/main'

// Thinking point: Do we want upd? meta? here

workflow PREPARE_GENS_INPUTS {
    take:
    ch_bam          // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_gvcf         // channel: [mandatory] [ val(meta), path(gvcf), path(gvcf_tbi), path(baf_positions) ]
    gatk_header     // channel: [mandatory] [ path(txt) ]

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
    MOSDEPTH(ch_mosdepth_in, ch_empty_fasta)

    // What output?
    // Some proper configuration needed here

    // Aha, need my own process here importing GAWK of course
    PREPROCESS_GENS_COV_INPUT(
        MOSDEPTH.out.regions_bed,
        gatk_header,
    )

    ch_gens_input = PREPROCESS_GENS_COV_INPUT.out.output.join(ch_gvcf)

    // Input: meta, coverage, gvcf
    GENERATE_GENS_DATA(
        ch_gens_input
    )

    // mosdepth
    // awk to GATK format
    // Python script

    // FIXME: Do we want a "without normal" cov output as well?
    emit:
    cov_bed_tbi = GENERATE_GENS_DATA.out.cov_bed_tbi
    baf_bed_tbi = GENERATE_GENS_DATA.out.baf_bed_tbi
    versions = ch_versions
}