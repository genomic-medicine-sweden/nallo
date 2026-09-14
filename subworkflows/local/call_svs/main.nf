include { DEBREAK_SV   } from '../debreak_sv/main'
include { HIFICNV_SV   } from '../hificnv_sv/main'
include { SAWFISH_SV   } from '../sawfish_sv/main'
include { SEVERUS_SV   } from '../severus_sv/main'
include { SNIFFLES1_SV } from '../sniffles1_sv/main'
include { SNIFFLES_SV  } from '../sniffles_sv/main'

workflow CALL_SVS {
    take:
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    ch_tandem_repeats // channel: [ val(meta), path(bed) ]
    ch_snvs // channel: [ val(meta), path(vcf) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_expected_xy_bed // channel: [ val(meta), path(bed) ]
    ch_expected_xx_bed // channel: [ val(meta), path(bed) ]
    ch_exclude_bed // channel: [ val(meta), path(bed) ]
    sv_callers_to_run //    List: [ 'caller1', 'caller2', ... ]
    force_sawfish_joint_call_single_samples //    bool: Force joint-calling with Sawfish even for single samples
    create_hificnv_maf_track //    bool: Should we create a MAF track for HiFiCNV calls?
    create_sawfish_maf_track //    bool: Should we create a MAF track for Sawfish calls?

    main:
    ch_sv_calls = channel.empty()

    if (sv_callers_to_run.contains('sniffles')) {
        SNIFFLES_SV(ch_bam_bai, ch_fasta, ch_tandem_repeats)
        ch_sv_calls = ch_sv_calls.mix(SNIFFLES_SV.out.vcf.map { meta, vcf -> [meta, vcf, []] })
    }

    if (sv_callers_to_run.contains('sniffles1')) {
        SNIFFLES1_SV(ch_bam_bai)
        ch_sv_calls = ch_sv_calls.mix(SNIFFLES1_SV.out.vcf.map { meta, vcf -> [meta, vcf, []] })
    }

    if (sv_callers_to_run.contains('severus')) {
        SEVERUS_SV(ch_bam_bai, ch_tandem_repeats)
        ch_sv_calls = ch_sv_calls.mix(SEVERUS_SV.out.vcf.map { meta, vcf -> [meta, vcf, []] })
    }

    if (sv_callers_to_run.contains('debreak')) {
        DEBREAK_SV(ch_bam_bai, ch_fasta)
        ch_sv_calls = ch_sv_calls.mix(DEBREAK_SV.out.vcf.map { meta, vcf -> [meta, vcf, []] })
    }

    if (sv_callers_to_run.contains('hificnv')) {
        HIFICNV_SV(ch_bam_bai, ch_snvs, ch_fasta, ch_expected_xy_bed, ch_expected_xx_bed, ch_exclude_bed, create_hificnv_maf_track)
        ch_sv_calls = ch_sv_calls.mix(
            HIFICNV_SV.out.vcf.join(HIFICNV_SV.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        )
    }

    if (sv_callers_to_run.contains('sawfish')) {
        SAWFISH_SV(ch_bam_bai, ch_snvs, ch_fasta, ch_expected_xy_bed, ch_expected_xx_bed, ch_exclude_bed, create_sawfish_maf_track, force_sawfish_joint_call_single_samples)
        ch_sv_calls = ch_sv_calls.mix(
            SAWFISH_SV.out.vcf.join(SAWFISH_SV.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        )
    }

    emit:
    sv_calls                           = ch_sv_calls // channel: [ val(meta), path(vcf), path(tbi) ] — tbi is [] for callers with skip_vep_prep: false
    hificnv_depth                      = sv_callers_to_run.contains('hificnv') ? HIFICNV_SV.out.depth : channel.empty() // channel: [ val(meta), path(bw) ]
    hificnv_copynum                    = sv_callers_to_run.contains('hificnv') ? HIFICNV_SV.out.copynum : channel.empty() // channel: [ val(meta), path(bedgraph) ]
    hificnv_maf                        = sv_callers_to_run.contains('hificnv') ? HIFICNV_SV.out.maf : channel.empty() // channel: [ val(meta), path(bw) ]
    sawfish_depth_bw                   = sv_callers_to_run.contains('sawfish') ? SAWFISH_SV.out.depth_bw : channel.empty() // channel: [ val(meta), path(bw) ]
    sawfish_copynum_bedgraph           = sv_callers_to_run.contains('sawfish') ? SAWFISH_SV.out.copynum_bedgraph : channel.empty() // channel: [ val(meta), path(bedgraph) ]
    sawfish_gc_bias_corrected_depth_bw = sv_callers_to_run.contains('sawfish') ? SAWFISH_SV.out.gc_bias_corrected_depth_bw : channel.empty() // channel: [ val(meta), path(bw) ]
    sawfish_maf_bw                     = sv_callers_to_run.contains('sawfish') ? SAWFISH_SV.out.maf_bw : channel.empty() // channel: [ val(meta), path(bw) ]
}
