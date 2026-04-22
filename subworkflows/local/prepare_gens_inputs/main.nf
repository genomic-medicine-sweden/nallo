include { CAT_CAT                       } from '../../../modules/nf-core/cat/cat/main'
include { GATK4_DENOISEREADCOUNTS       } from '../../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GAWK as MOSDEPTH_GATK_HEADER  } from '../../../modules/nf-core/gawk/main'
include { GAWK as MOSDEPTH_GATK_FORMAT  } from '../../../modules/nf-core/gawk/main'
include { MOSDEPTH                      } from '../../../modules/nf-core/mosdepth/main'
include { PREPARECOVANDBAF              } from '../../../modules/nf-core/gens/preparecovandbaf/main'
include { SAMTOOLS_VIEW                 } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT                 } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_AMPLICONCLIP         } from '../../../modules/nf-core/samtools/ampliconclip/main'
include { BAM_QC_FILTER                 } from '../../../modules/local/filter_bam/main'
include { GATK4_CLEANSAM               } from '../../../modules/local/gatk4/cleansam/main'

workflow PREPARE_GENS_INPUTS {
    take:
    ch_bam                     // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_gvcf                    // channel: [mandatory] [ val(meta), path(gvcfs)], [path(tbis) ]
    ch_baf_positions           // channel: [mandatory] [ val(meta), path(gz) ]
    ch_panel_of_normals_female // channel: [mandatory] [ val(meta), path(hd5) ]
    ch_panel_of_normals_male   // channel: [mandatory] [ val(meta), path(hd5) ]
    ch_mosdepth_bins           // channel: [mandatory] [ val(meta), path(bed) ]

    main:
/*    ch_bam
        .map { meta, bam, _bai -> [meta, bam] }
        .set { ch_bam_to_clip }

    // remove out of bounds reads (i.e. those that start within the chromosome in the reference but end beyond it) to avoid mosdepth errors
    SAMTOOLS_AMPLICONCLIP(
        ch_bam_to_clip,
        ch_mosdepth_bins
            .map { _meta, bed -> [bed] },
        false,
        false
    )

    SAMTOOLS_AMPLICONCLIP.out.bam.view()

    SAMTOOLS_SORT(
        SAMTOOLS_AMPLICONCLIP.out.bam,
        [[],[],[]],
        'bai'
    )

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_SORT.out.index)
        .set { ch_bam_bai_clipped }

    BAM_QC_FILTER(
        ch_bam,
    )

    BAM_QC_FILTER.out.bam_bai

    GATK4_CLEANSAM(
        ch_bam_to_clip,
        true
    )

    GATK4_CLEANSAM.out.bam_bai
        .set { ch_bam_bai_clipped }

    ch_bam_bai_clipped
        .combine(ch_mosdepth_bins)
        .map { meta, bam, bai, _bins_meta, bins ->
            [meta, bam, bai, bins]
        }
        .set { ch_mosdepth_in }

    ch_mosdepth_in.view()
*/

    ch_bam
        .combine(ch_mosdepth_bins)
        .map { meta, bam, bai, _bins_meta, bins ->
            [meta, bam, bai, bins]
        }
        .set { ch_mosdepth_in }


    // Prepare the header
    SAMTOOLS_VIEW(
        //ch_bam_bai_clipped,
        ch_bam,
        [[],[],[]],
        [],
        false
    )

    MOSDEPTH_GATK_HEADER(
        SAMTOOLS_VIEW.out.sam,
        [],
        false
    )

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
    MOSDEPTH_GATK_HEADER.out.output
        .join(MOSDEPTH_GATK_FORMAT.out.output)
        .map { meta, header, body -> [meta, [header, body]] }
        .set { ch_cat_input }

    CAT_CAT(ch_cat_input)

    CAT_CAT.out.file_out
        .branch { meta, _file ->
            male: meta.sex == 1
            female: meta.sex == 2
        }
        .set { ch_branched }

    ch_branched.male
        .combine(ch_panel_of_normals_male)
        .mix(ch_branched.female.combine(ch_panel_of_normals_female))
        .multiMap { meta, counts, _pon_meta, pon ->
            counts: [meta, counts]
            pon:    [meta, pon]
        }
        .set { ch_readcounts_input }


    // Calculate coverage
    GATK4_DENOISEREADCOUNTS(
        ch_readcounts_input.counts,
        ch_readcounts_input.pon,
    )

    GATK4_DENOISEREADCOUNTS.out.standardized
        .join(ch_gvcf)
        .set { ch_gens_input }

    // Generate final outputs
    ch_baf_positions
        .map { _meta, pos -> pos }
        .set { baf_positions }

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
}
