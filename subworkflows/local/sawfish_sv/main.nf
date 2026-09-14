include { SAWFISH_DISCOVER as RUN_SAWFISH_DISCOVER   } from '../../../modules/nf-core/sawfish/discover/main'
include { SAWFISH_JOINTCALL as RUN_SAWFISH_JOINTCALL } from '../../../modules/nf-core/sawfish/jointcall/main'

workflow SAWFISH_SV {
    take:
    ch_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_snvs // channel:  [optional] [ val(meta), path(vcf) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_expected_xy_bed // channel: [mandatory] [ val(meta), path(bed) ]
    ch_expected_xx_bed // channel: [mandatory] [ val(meta), path(bed) ]
    ch_exclude_bed // channel:  [optional] [ val(meta), path(bed) ]
    create_maf_track //     bool: [mandatory] should Sawfish produce a MAF track
    force_joint_call_single_samples //     bool: [mandatory] force joint-calling even for single samples

    main:
    ch_bam_vcf_for_discover = channel.empty()
    if (create_maf_track) {
        ch_bam_vcf_for_discover = ch_bam_bai.join(ch_snvs, failOnMismatch: true, failOnDuplicate: true)
    }
    else {
        ch_bam_vcf_for_discover = ch_bam_bai.map { meta, bam, bai -> [meta, bam, bai, []] }
    }

    ch_sawfish_discover_input = ch_bam_vcf_for_discover
        .combine(ch_expected_xx_bed)
        .combine(ch_expected_xy_bed)
        .multiMap { meta, bam, bai, vcf, xx_meta, xx_bed, xy_meta, xy_bed ->
            bam_bai: [meta, bam, bai]
            vcf: [meta, vcf]
            expected_copynumber_bed: meta.sex == 1
                ? [xy_meta, xy_bed]
                : [xx_meta, xx_bed]
        }

    RUN_SAWFISH_DISCOVER(
        ch_sawfish_discover_input.bam_bai,
        ch_fasta,
        ch_sawfish_discover_input.expected_copynumber_bed,
        ch_sawfish_discover_input.vcf,
        ch_exclude_bed,
    )

    ch_sawfish_jointcall_input = RUN_SAWFISH_DISCOVER.out.discover_dir
        .join(ch_sawfish_discover_input.bam_bai, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, discover_dir, bam, bai ->
            def new_meta = force_joint_call_single_samples
                ? meta
                : [id: meta.family_id, family_id: meta.family_id]
            [new_meta, discover_dir, bam, bai]
        }
        .groupTuple()
        .multiMap { meta, discover_dirs, bams, bais ->
            dir: [meta, discover_dirs]
            bam_bai: [meta, bams, bais]
        }

    RUN_SAWFISH_JOINTCALL(
        ch_sawfish_jointcall_input.dir,
        ch_fasta,
        ch_sawfish_jointcall_input.bam_bai,
        [[], []],
    )

    ch_vcf_tbi = RUN_SAWFISH_JOINTCALL.out.vcf
        .join(RUN_SAWFISH_JOINTCALL.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi -> [meta + [sv_caller: 'sawfish', needs_reheader: false, skip_vep_prep: true], vcf, tbi] }

    emit:
    vcf                        = ch_vcf_tbi.map { meta, vcf, _tbi -> [meta, vcf] } // channel: [ val(meta), path(vcf) ]
    tbi                        = ch_vcf_tbi.map { meta, _vcf, tbi -> [meta, tbi] } // channel: [ val(meta), path(tbi) ]
    depth_bw                   = addSampleIdFromSawfishPath(RUN_SAWFISH_JOINTCALL.out.depth_bw) // channel: [ val(meta), path(bw) ]
    copynum_bedgraph           = addSampleIdFromSawfishPath(RUN_SAWFISH_JOINTCALL.out.copynum_bedgraph) // channel: [ val(meta), path(bedgraph) ]
    gc_bias_corrected_depth_bw = addSampleIdFromSawfishPath(RUN_SAWFISH_JOINTCALL.out.gc_bias_corrected_depth_bw) // channel: [ val(meta), path(bw) ]
    maf_bw                     = addSampleIdFromSawfishPath(RUN_SAWFISH_JOINTCALL.out.maf_bw) // channel: [ val(meta), path(bw) ]
}

def addSampleIdFromSawfishPath(ch_sawfish_output) {
    ch_sawfish_output
        .transpose()
        .map { meta, file ->
            def sample_id = file.parent.name.replaceFirst(/[^_]*_/, "")
            [meta + [id: sample_id], file]
        }
}
