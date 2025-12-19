include { BCFTOOLS_CONCAT }               from '../../../modules/nf-core/bcftools/concat/main'
include { MOSDEPTH }                      from '../../../modules/nf-core/mosdepth/main'
include { GENERATE_GENS_DATA }            from '../../../modules/local/generate_gens_data/main'
include { MOSDEPTH_TO_GATK_FORMAT }       from '../../../modules/local/mosdepth_to_gatk_format/main'
include { GATK4_DENOISEREADCOUNTS }       from '../../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GAWK }                          from '../../../modules/nf-core/gawk/main'
include { SAMTOOLS_VIEW }                 from '../../../modules/nf-core/samtools/view/main'

workflow PREPARE_GENS_INPUTS {
    take:
    ch_bam              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_gvcfs            // channel: [mandatory] [ val(meta), [path(gvcfs)], [path(tbis)] ]
    baf_positions       // value:   [mandatory] [ path(gz) ]
    panel_of_normals    // value:   [optional]  [ path(hd5) ]

    main:
    ch_versions = channel.empty()

    BCFTOOLS_CONCAT(
        ch_gvcfs
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    BCFTOOLS_CONCAT.out.vcf
        .join(BCFTOOLS_CONCAT.out.tbi)
        .set { ch_vcf_tbi }

    ch_bam
        .map { meta, bam, bai -> tuple(meta, bam, bai, []) }
        .set { ch_mosdepth_in }

    MOSDEPTH(
        ch_mosdepth_in,
        [[],[]]
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    SAMTOOLS_VIEW(
        ch_bam,
        [[],[]],
        [],
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    GAWK(
        SAMTOOLS_VIEW.out.sam,
        [],
        false
    )
    ch_versions = ch_versions.mix(GAWK.out.versions)

    MOSDEPTH.out.regions_bed
        .join(GAWK.out.output)
        .set { ch_mosdepth_to_gatk_in }

    MOSDEPTH_TO_GATK_FORMAT(
        ch_mosdepth_to_gatk_in
    )
    ch_versions = ch_versions.mix(MOSDEPTH_TO_GATK_FORMAT.out.versions)

    MOSDEPTH_TO_GATK_FORMAT.out.output
        .map { meta, _tsv -> [ meta, panel_of_normals ] }
        .set { ch_pon }

    GATK4_DENOISEREADCOUNTS(
        MOSDEPTH_TO_GATK_FORMAT.out.output,
        ch_pon
    )
    ch_versions = ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions)

        GATK4_DENOISEREADCOUNTS.out.standardized
        .set { ch_cov }

    ch_cov
        .join(ch_vcf_tbi)
        .set { ch_gens_input }

    GENERATE_GENS_DATA(
        ch_gens_input,
        baf_positions
    )
    ch_versions = ch_versions.mix(GENERATE_GENS_DATA.out.versions)

    emit:
    cov_bed_tbi = GENERATE_GENS_DATA.out.cov_bed_tbi // channel: [ val(meta), path(bed_gz), path(tbi) ]
    baf_bed_tbi = GENERATE_GENS_DATA.out.baf_bed_tbi // channel: [ val(meta), path(bed_gz), path(tbi) ]
    versions = ch_versions                           // channel: [ path(versions.yml) ]
}
