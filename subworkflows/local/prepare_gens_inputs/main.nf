include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'
include { GENERATE_GENS_DATA } from '../../../modules/local/generate_gens_data/main'
include { MOSDEPTH_TO_GATK_FORMAT } from '../../../modules/local/mosdepth_to_gatk_format/main'
include { GATK4_DENOISEREADCOUNTS } from '../../../modules/nf-core/gatk4/denoisereadcounts/main'

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

    MOSDEPTH(ch_mosdepth_in, ch_empty_fasta)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    MOSDEPTH_TO_GATK_FORMAT(
        MOSDEPTH.out.regions_bed,
        gatk_header,
    )
    ch_versions = ch_versions.mix(MOSDEPTH_TO_GATK_FORMAT.out.versions)

    // BGZIP
    // Tabix

    // What output?
    // Some proper configuration needed here

    channel
        .fromPath(panel_of_normals)
        .map { f -> tuple([], f) }
        .set { ch_pon }

    GATK4_DENOISEREADCOUNTS(
        MOSDEPTH_TO_GATK_FORMAT.out.output,
        ch_pon
        // meta, pon
    )
    ch_versions = ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions)

    // PlotDenoisedCopyRatios
    // Cleanup

    // Consider using GATK4_DENOISEREADCOUNTS.out.standardized as well here
    // Also worth considering showing the raw data or normalized raw
    // without the PoN

    // FIXME: What output
    ch_gens_input = GATK4_DENOISEREADCOUNTS.out.standardized.join(ch_gvcf)

    ch_gens_input.view()

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
