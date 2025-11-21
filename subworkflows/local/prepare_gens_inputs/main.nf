include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'
include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'
include { GENERATE_GENS_DATA } from '../../../modules/local/generate_gens_data/main'
include { MOSDEPTH_TO_GATK_FORMAT } from '../../../modules/local/mosdepth_to_gatk_format/main'
include { GATK4_DENOISEREADCOUNTS } from '../../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GENERATE_MOSDEPTH_GATK_HEADER } from '../../../modules/local/generate_mosdepth_gatk_header/main.nf'
include { NORMALIZE_MOSDEPTH_COVERAGE } from '../../../modules/local/normalize_mosdepth_coverage/main.nf'

workflow PREPARE_GENS_INPUTS {
    take:
    ch_bam              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_gvcfs            // channel: [mandatory] [ val(meta), tuple(path(gvcfs)), tuple(path(tbis))]
    baf_positions       // value:   [mandatory] [ path(gz) ]
    panel_of_normals    // value:   [optional] [ path(hd5) ]
    use_pon             // value:   [mandatory] [ boolean ]

    main:
    ch_versions = channel.empty()

    BCFTOOLS_CONCAT(ch_gvcfs)
    ch_vcf_tbi = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi)

    ch_mosdepth_in = ch_bam.map { meta, bam, bai -> tuple(meta, bam, bai, []) }
    ch_empty_fasta = ch_bam.map { meta, _bam, _bai -> tuple(meta, []) }

    MOSDEPTH(ch_mosdepth_in, ch_empty_fasta)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    GENERATE_MOSDEPTH_GATK_HEADER(ch_bam)
    ch_versions = ch_versions.mix(GENERATE_MOSDEPTH_GATK_HEADER.out.versions)
    MOSDEPTH.out.regions_bed
        .join(GENERATE_MOSDEPTH_GATK_HEADER.out.header)
        .set { ch_mosdepth_to_gatk_in }
    MOSDEPTH_TO_GATK_FORMAT(ch_mosdepth_to_gatk_in)
    ch_versions = ch_versions.mix(MOSDEPTH_TO_GATK_FORMAT.out.versions)

    // PON path
    if (use_pon) {
        MOSDEPTH_TO_GATK_FORMAT.out.output
            .map { meta, _tsv -> tuple(meta, panel_of_normals) }
            .set { ch_pon }
        GATK4_DENOISEREADCOUNTS(
            MOSDEPTH_TO_GATK_FORMAT.out.output,
            ch_pon
        )
        ch_versions = ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions)
        GATK4_DENOISEREADCOUNTS.out.standardized
            .set { ch_cov }
    } else {
        NORMALIZE_MOSDEPTH_COVERAGE(MOSDEPTH_TO_GATK_FORMAT.out.output)
        NORMALIZE_MOSDEPTH_COVERAGE.out.normalized
            .set { ch_cov }
    }

    ch_gens_input = ch_cov.join(ch_vcf_tbi)

    GENERATE_GENS_DATA(
        ch_gens_input,
        baf_positions
    )
    ch_versions = ch_versions.mix(GENERATE_GENS_DATA.out.versions)

    emit:
    cov_bed_tbi = GENERATE_GENS_DATA.out.cov_bed_tbi
    baf_bed_tbi = GENERATE_GENS_DATA.out.baf_bed_tbi
    versions = ch_versions
}
