include { CRAMINO        } from '../../../modules/nf-core/cramino/main'
include { FASTQC         } from '../../../modules/nf-core/fastqc/main'
include { MOSDEPTH       } from '../../../modules/nf-core/mosdepth/main'
include { SAMBAMBA_DEPTH } from '../../../modules/nf-core/sambamba/depth/main'
workflow QC_ALIGNED_READS {
    take:
    ch_bam_bai         // channel: [ val(meta), [bam, bai] ]
    ch_fasta           // channel: [ val(meta), fasta ]
    ch_mosdepth_bed    // channel: [ val(meta), bed ]
    ch_sambamba_bed    // channel: [ val(meta), bed ]
    run_sambamba_depth //    bool: Should sambamba depth be run?

    main:
    ch_sambamba_depth_bed = channel.empty()
    ch_versions           = channel.empty()

    FASTQC(
        ch_bam_bai.map { meta, bam, _bai -> [meta, bam] }
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    CRAMINO(
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(CRAMINO.out.versions)

    ch_bam_bai
        .combine(ch_mosdepth_bed.map { _meta, bed -> bed }.toList()) // toList() enables passing [] if ch_bed is empty
        .set { mosdepth_in }

    MOSDEPTH(
        mosdepth_in,
        ch_fasta,
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    if (run_sambamba_depth) {
        SAMBAMBA_DEPTH(
            ch_bam_bai,
            ch_sambamba_bed,
            'region',
        )
        ch_versions = ch_versions.mix(SAMBAMBA_DEPTH.out.versions)

        SAMBAMBA_DEPTH.out.bed
            .set { ch_sambamba_depth_bed }
    }

    emit:
    fastqc_zip           = FASTQC.out.zip           // channel: [ val(meta), path(zip) ]
    mosdepth_summary     = MOSDEPTH.out.summary_txt // channel: [ val(meta), path(txt) ]
    mosdepth_global_dist = MOSDEPTH.out.global_txt  // channel: [ val(meta), path(txt) ]
    mosdepth_region_dist = MOSDEPTH.out.regions_txt // channel: [ val(meta), path(txt) ]
    sambamba_depth_bed   = ch_sambamba_depth_bed    // channel: [ val(meta), path(bed) ]
    versions             = ch_versions              // channel: [ versions.yml ]
}
