include { CRAMINO        } from '../../../modules/nf-core/cramino/main'
include { FASTQC         } from '../../../modules/nf-core/fastqc/main'
include { MOSDEPTH       } from '../../../modules/nf-core/mosdepth/main'
include { SAMBAMBA_DEPTH } from '../../../modules/nf-core/sambamba/depth/main'
include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view/main'
workflow QC_ALIGNED_READS {
    take:
    ch_bam_bai // channel: [ val(meta), [bam, bai] ]
    ch_fasta // channel: [ val(meta), fasta ]
    ch_mosdepth_bed // channel: [ val(meta), bed ]
    ch_sambamba_bed // channel: [ val(meta), bed ]
    run_sambamba_depth // bool: Should sambamba depth be run?
    ch_cramino_bed // channel: [ val(meta), bed ] or empty; if supplied, reads are filtered to regions in the bed file before cramino QC

    main:
    ch_sambamba_depth_bed = channel.empty()

    FASTQC(
        ch_bam_bai.map { meta, bam, _bai -> [meta, bam] }
    )

    ch_bam_bai
        .combine(ch_cramino_bed.map { _meta, bed -> bed }.toList())
        .branch { meta, bam, bai, bed ->
            bed: bed
            no_bed: !bed
        }
        .set { ch_cramino_branches }

    SAMTOOLS_VIEW(
        ch_cramino_branches.bed.map { meta, bam, bai, _bed -> [meta, bam, bai] },
        [[], [], []],
        [[], []],
        ch_cramino_branches.bed.map { meta, _bam, _bai, bed -> [meta, bed] },
        'bai',
    )

    ch_bam_bai_for_cramino = ch_cramino_branches.no_bed
        .map { meta, bam, bai, _bed -> [meta, bam, bai] }
        .mix(
            SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_VIEW.out.bai)
        )

    CRAMINO(
        ch_bam_bai_for_cramino
    )

    // toList() enables passing [] if ch_bed is empty
    ch_mosdepth_in = ch_bam_bai.combine(ch_mosdepth_bed.map { _meta, bed -> bed }.toList())

    MOSDEPTH(
        ch_mosdepth_in,
        ch_fasta,
    )

    if (run_sambamba_depth) {
        SAMBAMBA_DEPTH(
            ch_bam_bai,
            ch_sambamba_bed,
            'region',
        )

        ch_sambamba_depth_bed = SAMBAMBA_DEPTH.out.bed
    }

    emit:
    cramino_stats           = CRAMINO.out.stats // channel: [ val(meta), path(txt)        ]
    cramino_arrow           = CRAMINO.out.arrow // channel: [ val(meta), path(arrow)      ]
    fastqc_html             = FASTQC.out.html // channel: [ val(meta), path(html)       ]
    fastqc_zip              = FASTQC.out.zip // channel: [ val(meta), path(zip)        ]
    mosdepth_summary        = MOSDEPTH.out.summary_txt // channel: [ val(meta), path(txt)        ]
    mosdepth_global_dist    = MOSDEPTH.out.global_txt // channel: [ val(meta), path(txt)        ]
    mosdepth_regions_dist   = MOSDEPTH.out.regions_txt // channel: [ val(meta), path(txt)        ]
    mosdepth_per_base_d4    = MOSDEPTH.out.per_base_d4 // channel: [ val(meta), path(d4)         ]
    mosdepth_regions_bed    = MOSDEPTH.out.regions_bed // channel: [ val(meta), path(bed.gz)     ]
    mosdepth_regions_csi    = MOSDEPTH.out.regions_csi // channel: [ val(meta), path(bed.gz.csi) ]
    mosdepth_thresholds_bed = MOSDEPTH.out.thresholds_bed // channel: [ val(meta), path(bed.gz)     ]
    mosdepth_thresholds_csi = MOSDEPTH.out.thresholds_csi // channel: [ val(meta), path(bed.gz.csi) ]
    sambamba_depth_bed      = ch_sambamba_depth_bed // channel: [ val(meta), path(bed)        ]
}
