include { BCFTOOLS_CONCAT }               from '../../../modules/nf-core/bcftools/concat/main'
include { CAT_CAT }                       from '../../../modules/nf-core/cat/cat/main'
include { GATK4_DENOISEREADCOUNTS }       from '../../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GAWK as MOSDEPTH_GATK_HEADER }  from '../../../modules/nf-core/gawk/main'
include { GAWK as MOSDEPTH_GATK_FORMAT }  from '../../../modules/nf-core/gawk/main'
include { MOSDEPTH }                      from '../../../modules/nf-core/mosdepth/main'
include { PREPARECOVANDBAF }              from '../../../modules/nf-core/gens/preparecovandbaf/main'
include { SAMTOOLS_VIEW }                 from '../../../modules/nf-core/samtools/view/main'
include { TABIX_BGZIP }                   from '../../../modules/nf-core/tabix/bgzip/main'

workflow PREPARE_GENS_INPUTS {
    take:
    ch_bam              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_gvcf             // channel: [mandatory] [ val(meta), path(gvcfs)], [path(tbis) ]
    baf_positions       // value:   [mandatory] [ path(gz) ]
    panel_of_normals    // value:   [mandatory] [ path(hd5) ]
    mosdepth_bins       // value:   [mandatory] [ path(bed) ]

    main:
    ch_versions = channel.empty()

    ch_bam
        .map { meta, bam, bai -> tuple(meta, bam, bai, mosdepth_bins) }
        .set { ch_mosdepth_in }

    // Prepare the header
    SAMTOOLS_VIEW(
        ch_bam,
        [[],[]],
        [],
        false
    )
    ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    MOSDEPTH_GATK_HEADER(
        SAMTOOLS_VIEW.out.sam,
        [],
        false
    )
    ch_versions.mix(MOSDEPTH_GATK_HEADER.out.versions)

    // Prepare the body
    MOSDEPTH(
        ch_mosdepth_in,
        [[],[]]
    )
    ch_versions.mix(MOSDEPTH.out.versions)
    MOSDEPTH_GATK_FORMAT(
        MOSDEPTH.out.regions_bed,
        [],
        false
    )
    ch_versions.mix(MOSDEPTH_GATK_FORMAT.out.versions)

    // Prepare GATK inputs
    MOSDEPTH_GATK_HEADER.out.output
        .join(MOSDEPTH_GATK_FORMAT.out.output)
        .map { meta, header, body -> tuple(meta, [header, body]) }
        .set { ch_cat_input }
    CAT_CAT(ch_cat_input)
    ch_versions.mix(CAT_CAT.out.versions)
    CAT_CAT.out.file_out
        .map { meta, _tsv -> [ meta, panel_of_normals ] }
        .set { ch_pon }

    // Calculate coverage
    GATK4_DENOISEREADCOUNTS(
        CAT_CAT.out.file_out,
        ch_pon
    )
    ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions)
    GATK4_DENOISEREADCOUNTS.out.standardized
        .set { ch_cov }
    ch_cov
        .join(ch_gvcf)
        .set { ch_gens_input }

    // Generate final outputs
    PREPARECOVANDBAF(
        ch_gens_input,
        baf_positions
    )

    ch_cov_gz_tbi = PREPARECOVANDBAF.out.cov_gz
        .join(PREPARECOVANDBAF.out.cov_tbi)
    ch_baf_gz_tbi = PREPARECOVANDBAF.out.baf_gz
        .join(PREPARECOVANDBAF.out.baf_tbi)

    emit:
    cov_bed_tbi = ch_cov_gz_tbi    // channel: [ val(meta), path(bed_gz), path(tbi) ]
    baf_bed_tbi = ch_baf_gz_tbi    // channel: [ val(meta), path(bed_gz), path(tbi) ]
    versions = ch_versions         // channel: [ path(versions.yml) ]
}
