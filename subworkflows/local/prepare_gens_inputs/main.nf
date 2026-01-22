include { BCFTOOLS_CONCAT }               from '../../../modules/nf-core/bcftools/concat/main'
include { MOSDEPTH }                      from '../../../modules/nf-core/mosdepth/main'
include { PREPARECOVANDBAF }              from '../../../modules/nf-core/gens/preparecovandbaf/main'
include { MOSDEPTH_TO_GATK_FORMAT }       from '../../../modules/local/mosdepth_to_gatk_format/main'
include { GATK4_DENOISEREADCOUNTS }       from '../../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GAWK }                          from '../../../modules/nf-core/gawk/main'
include { SAMTOOLS_VIEW }                 from '../../../modules/nf-core/samtools/view/main'

workflow PREPARE_GENS_INPUTS {
    take:
    ch_bam              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_gvcf             // channel: [mandatory] [ val(meta), path(gvcfs)], [path(tbis) ]
    baf_positions       // value:   [mandatory] [ path(gz) ]
    panel_of_normals    // value:   [optional]  [ path(hd5) ]

    main:

    ch_versions = channel.empty()

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

    // FIXME: Replace this one with built-in commands?
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
        .join(ch_gvcf)
        .set { ch_gens_input }

    print(">>> Before prepare gens input data")

    PREPARECOVANDBAF(
        ch_gens_input,
        baf_positions
    )

    print(">>> After prepare gens input data")
    PREPARECOVANDBAF.out.cov_gz.view()
    PREPARECOVANDBAF.out.cov_tbi.view()

    ch_versions = ch_versions.mix(PREPARECOVANDBAF.out.versions)

    ch_cov_gz_tbi = PREPARECOVANDBAF.out.cov_gz.join(PREPARECOVANDBAF.out.cov_tbi)
    ch_baf_gz_tbi = PREPARECOVANDBAF.out.baf_gz.join(PREPARECOVANDBAF.out.baf_tbi)

    ch_cov_gz_tbi.view()
    ch_baf_gz_tbi.view()

    emit:
    cov_bed_tbi = ch_cov_gz_tbi                      // channel: [ val(meta), path(bed_gz), path(tbi) ]
    baf_bed_tbi = ch_baf_gz_tbi                      // channel: [ val(meta), path(bed_gz), path(tbi) ]
    versions = ch_versions                           // channel: [ path(versions.yml) ]
}
