include { BCFTOOLS_CONCAT }               from '../../../modules/nf-core/bcftools/concat/main'
include { MOSDEPTH }                      from '../../../modules/nf-core/mosdepth/main'
include { PREPARECOVANDBAF }              from '../../../modules/nf-core/gens/preparecovandbaf/main'
include { GATK4_DENOISEREADCOUNTS }       from '../../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GAWK as MOSDEPTH_GATK_HEADER }  from '../../../modules/nf-core/gawk/main'
include { GAWK as MOSDEPTH_GATK_FORMAT }  from '../../../modules/nf-core/gawk/main'
include { SAMTOOLS_VIEW }                 from '../../../modules/nf-core/samtools/view/main'
include { CAT_CAT }                       from '../../../modules/nf-core/cat/cat/main'
include { TABIX_BGZIP }                   from '../../../modules/nf-core/tabix/bgzip/main'

workflow PREPARE_GENS_INPUTS {
    take:
    ch_bam              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_gvcf             // channel: [mandatory] [ val(meta), path(gvcfs)], [path(tbis) ]
    baf_positions       // value:   [mandatory] [ path(gz) ]
    panel_of_normals    // value:   [optional]  [ path(hd5) ]

    main:

    ch_bam
        .map { meta, bam, bai -> tuple(meta, bam, bai, []) }
        .set { ch_mosdepth_in }

    // Prepare the header
    SAMTOOLS_VIEW(
        ch_bam,
        [[],[]],
        [],
        false
    )
    MOSDEPTH_GATK_HEADER(
        SAMTOOLS_VIEW.out.sam,
        [],
        false
    )
    TABIX_BGZIP(MOSDEPTH_GATK_HEADER.out.output)

    // Prepare the body
    MOSDEPTH(
        ch_mosdepth_in,
        [[],[]]
    )
    MOSDEPTH_GATK_FORMAT(
        MOSDEPTH.out.regions_bed,
        [],
        false
    )

    // Prepare GATK inputs
    TABIX_BGZIP.out.output
        .join(MOSDEPTH_GATK_FORMAT.out.output)
        .map { meta, header, body -> tuple(meta, [header, body]) }
        .set { ch_cat_input }
    CAT_CAT(ch_cat_input)
    CAT_CAT.out.file_out
        .map { meta, _tsv -> [ meta, panel_of_normals ] }
        .set { ch_pon }

    // Calculate coverage
    GATK4_DENOISEREADCOUNTS(
        CAT_CAT.out.file_out,
        ch_pon
    )
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

    ch_cov_gz_tbi = PREPARECOVANDBAF.out.cov_gz.join(PREPARECOVANDBAF.out.cov_tbi)
    ch_baf_gz_tbi = PREPARECOVANDBAF.out.baf_gz.join(PREPARECOVANDBAF.out.baf_tbi)

    emit:
    cov_bed_tbi = ch_cov_gz_tbi                      // channel: [ val(meta), path(bed_gz), path(tbi) ]
    baf_bed_tbi = ch_baf_gz_tbi                      // channel: [ val(meta), path(bed_gz), path(tbi) ]
}
