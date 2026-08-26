include { BCFTOOLS_VIEW                     } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_SORT                     } from '../../../modules/nf-core/bcftools/sort/main'
include { DEBREAK                           } from '../../../modules/nf-core/debreak/main'
include { HIFICNV                           } from '../../../modules/nf-core/hificnv/main'
include { SAWFISH_DISCOVER                  } from '../../../modules/nf-core/sawfish/discover/main'
include { SAWFISH_JOINTCALL                 } from '../../../modules/nf-core/sawfish/jointcall/main'
include { SEVERUS                           } from '../../../modules/nf-core/severus/main'
include { SNIFFLES                          } from '../../../modules/nf-core/sniffles/main'
include { SNIFFLES1                         } from '../../../modules/local/sniffles1/main'
include { TABIX_TABIX as TABIX_HIFICNV      } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as TABIX_SEVERUS } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { VEP_PREP_SV                       } from '../../../modules/local/vep_prep_sv/main'

workflow CALL_SVS {
    take:
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    ch_tandem_repeats // channel: [ val(meta), path(bed) ]
    ch_snvs // channel: [ val(meta), path(vcf) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_expected_xy_bed // channel: [ val(meta), path(bed) ]
    ch_expected_xx_bed // channel: [ val(meta), path(bed) ]
    ch_exclude_bed // channel: [ val(meta), path(bed) ]
    sv_callers_to_run //    List: [ 'caller1', 'caller2', 'caller3' ]
    ch_sv_call_regions // channel: [ val(meta), path(bed) ]
    filter_calls_on_regions //    bool: Should we filter SV calls to the regions provided in ch_sv_call_regions?
    force_sawfish_joint_call_single_samples //    bool: Force joint-calling with Sawfish even for single samples
    create_hificnv_maf_track //    bool: Should we create a MAF track for HiFiCNV/Sawfish calls?
    create_sawfish_maf_track //    bool: Should we create a MAF track for HiFiCNV/Sawfish calls?

    main:
    ch_sv_calls = channel.empty()
    ch_for_vep_prep_sv = channel.empty()

    //
    // Call SVs with Severus
    //
    if (sv_callers_to_run.contains('severus')) {

        SEVERUS(
            ch_bam_bai.map { meta, bam, bai -> [meta, bam, bai, [], [], []] },
            ch_tandem_repeats,
        )

        ch_for_vep_prep_sv = ch_for_vep_prep_sv.mix(
            SEVERUS.out.all_vcf.map { meta, vcf -> [meta + [sv_caller: 'severus'], vcf] }
        )
    }

    //
    // Call SVs with Sniffles
    //
    if (sv_callers_to_run.contains('sniffles')) {

        SNIFFLES(
            ch_bam_bai,
            ch_fasta,
            ch_tandem_repeats,
            true,
            false,
        )

        ch_for_vep_prep_sv = ch_for_vep_prep_sv.mix(
            SNIFFLES.out.vcf.map { meta, vcf -> [meta + [sv_caller: 'sniffles'], vcf] }
        )
    }

    //
    // Call SVs with Sniffles v1
    //
    if (sv_callers_to_run.contains('sniffles1')) {

        SNIFFLES1(
            ch_bam_bai
        )

        ch_for_vep_prep_sv = ch_for_vep_prep_sv.mix(
            SNIFFLES1.out.vcf.map { meta, vcf -> [meta + [sv_caller: 'sniffles1'], vcf] }
        )
    }

    //
    // Call SVs with DeBreak
    //
    if (sv_callers_to_run.contains('debreak')) {

        DEBREAK(
            ch_bam_bai,
            ch_fasta,
        )

        ch_for_vep_prep_sv = ch_for_vep_prep_sv.mix(
            DEBREAK.out.vcf.map { meta, vcf -> [meta + [sv_caller: 'debreak'], vcf] }
        )
    }

    //
    // Prepare SV VCFs for annotation
    //
    VEP_PREP_SV(ch_for_vep_prep_sv)

    if (sv_callers_to_run.contains('severus')) {

        TABIX_SEVERUS(
            VEP_PREP_SV.out.vcf.filter { meta, _vcf -> meta.sv_caller == 'severus' }
        )

        ch_sv_calls = ch_sv_calls.mix(TABIX_SEVERUS.out.gz_index)
    }

    BCFTOOLS_SORT(
        VEP_PREP_SV.out.vcf.filter { meta, _vcf -> meta.sv_caller in ['sniffles', 'sniffles1', 'debreak'] }
    )

    ch_sv_calls = ch_sv_calls.mix(
        BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_SORT.out.tbi, failOnMismatch: true, failOnDuplicate: true)
    )

    //
    // Call CNVs with HiFiCNV
    //
    if (sv_callers_to_run.contains('hificnv')) {

        // Join SNV VCFs into input channels only if we want MAF track
        // Otherwise, we can skip the join and just pass an empty list to the module, since the MAF track is the only thing that needs the SNV VCFs.
        if (create_hificnv_maf_track) {
            ch_for_hificnv = ch_bam_bai.join(ch_snvs, failOnMismatch: true, failOnDuplicate: true)
        }
        else {
            ch_for_hificnv = ch_bam_bai.map { meta, bam, bai -> [meta, bam, bai, []] }
        }

        // Select expected copynumber BED based on sex before passing it to the module
        ch_hificnv_input = ch_for_hificnv
            .combine(ch_expected_xy_bed)
            .combine(ch_expected_xx_bed)
            .multiMap { meta, bam, bai, maf, xy_meta, xy_bed, xx_meta, xx_bed ->
                def expected_cn_meta = meta.sex == 1 ? xy_meta : xx_meta
                def expected_cn_bed = meta.sex == 1 ? xy_bed : xx_bed
                bam_bai_maf: [meta, bam, bai, maf]
                expected_cn: [expected_cn_meta, expected_cn_bed]
            }

        HIFICNV(
            ch_hificnv_input.bam_bai_maf,
            ch_fasta,
            ch_exclude_bed,
            ch_hificnv_input.expected_cn,
        )

        TABIX_HIFICNV(
            HIFICNV.out.vcf
        )

        ch_sv_calls = ch_sv_calls.mix(
            addCallerToMeta(
                HIFICNV.out.vcf.join(TABIX_HIFICNV.out.index, failOnMismatch: true, failOnDuplicate: true),
                'hificnv',
            )
        )
    }

    //
    // Call SVs with Sawfish
    //
    if (sv_callers_to_run.contains('sawfish')) {

        // Join SNV VCFs into input channels only if we want MAF track
        // Otherwise, we can skip the join and just pass an empty list to the module, since the MAF track is the only thing that needs the SNV VCFs.
        if (create_sawfish_maf_track) {
            ch_bam_vcf_for_sawfish_discover = ch_bam_bai.join(ch_snvs, failOnMismatch: true, failOnDuplicate: true)
        }
        else {
            ch_bam_vcf_for_sawfish_discover = ch_bam_bai.map { meta, bam, bai -> [meta, bam, bai, []] }
        }

        ch_sawfish_discover_input = ch_bam_vcf_for_sawfish_discover
            .combine(ch_expected_xx_bed)
            .combine(ch_expected_xy_bed)
            .multiMap { meta, bam, bai, vcf, xx_meta, xx_bed, xy_meta, xy_bed ->
                bam_bai: [meta, bam, bai]
                vcf: [meta, vcf]
                expected_copynumber_bed: meta.sex == 1
                    ? [xy_meta, xy_bed]
                    : [xx_meta, xx_bed]
            }

        SAWFISH_DISCOVER(
            ch_sawfish_discover_input.bam_bai,
            ch_fasta,
            ch_sawfish_discover_input.expected_copynumber_bed,
            ch_sawfish_discover_input.vcf,
            ch_exclude_bed,
        )

        // Sawfish needs joint-calling to actually produce SV calls. Without it, there are no sample names
        // in the VCFs, and they can't be post-processed with bcftools. Therefore, we do joint-calling step
        // here directly, and skip doing it later with SVDB merging.
        ch_sawfish_jointcall_input = SAWFISH_DISCOVER.out.discover_dir
            .join(ch_sawfish_discover_input.bam_bai, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, discover_dir, bam, bai ->
                def new_meta = force_sawfish_joint_call_single_samples
                    ? meta
                    : [id: meta.family_id, family_id: meta.family_id]

                [new_meta, discover_dir, bam, bai]
            }
            .groupTuple()
            .multiMap { meta, discover_dirs, bams, bais ->
                dir: [meta, discover_dirs]
                bam_bai: [meta, bams, bais]
            }

        SAWFISH_JOINTCALL(
            ch_sawfish_jointcall_input.dir,
            ch_fasta,
            ch_sawfish_jointcall_input.bam_bai,
            [[], []],
        )

        ch_sv_calls = ch_sv_calls.mix(
            addCallerToMeta(
                SAWFISH_JOINTCALL.out.vcf.join(SAWFISH_JOINTCALL.out.tbi, failOnMismatch: true, failOnDuplicate: true),
                'sawfish',
            )
        )
    }

    //
    // Post-process SV calls
    //
    if (filter_calls_on_regions) {

        BCFTOOLS_VIEW(
            ch_sv_calls,
            ch_sv_call_regions.map { _meta, bed -> bed },
            [],
            [],
        )

        ch_sv_calls_filtered = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi, failOnMismatch: true, failOnDuplicate: true)
    }
    else {
        ch_sv_calls_filtered = ch_sv_calls
    }

    emit:
    vcf                                = ch_sv_calls_filtered // channel: [ val(meta),  path(vcf), path(tbi) ]
    hificnv_depth                      = sv_callers_to_run.contains('hificnv') ? HIFICNV.out.depth : channel.empty() // channel: [ val(meta), path(bw) ]
    hificnv_copynum                    = sv_callers_to_run.contains('hificnv') ? HIFICNV.out.copynum : channel.empty() // channel: [ val(meta), path(bedgraph) ]
    hificnv_maf                        = sv_callers_to_run.contains('hificnv') ? HIFICNV.out.maf : channel.empty() // channel: [ val(meta), path(bw) ]
    sawfish_depth_bw                   = sv_callers_to_run.contains('sawfish') ? addSampleIdFromSawfishPath(SAWFISH_JOINTCALL.out.depth_bw) : channel.empty() // channel: [ val(meta), path(bw) ]
    sawfish_copynum_bedgraph           = sv_callers_to_run.contains('sawfish') ? addSampleIdFromSawfishPath(SAWFISH_JOINTCALL.out.copynum_bedgraph) : channel.empty() // channel: [ val(meta), path(bedgraph) ]
    sawfish_gc_bias_corrected_depth_bw = sv_callers_to_run.contains('sawfish') ? addSampleIdFromSawfishPath(SAWFISH_JOINTCALL.out.gc_bias_corrected_depth_bw) : channel.empty() // channel: [ val(meta), path(bw) ]
    sawfish_maf_bw                     = sv_callers_to_run.contains('sawfish') ? addSampleIdFromSawfishPath(SAWFISH_JOINTCALL.out.maf_bw) : channel.empty() // channel: [ val(meta), path(bw) ]
}

def addCallerToMeta(ch_caller_calls, sv_caller) {
    ch_caller_calls.map { meta, vcf, tbi ->
        [meta + [sv_caller: sv_caller], vcf, tbi]
    }
}

def addSampleIdFromSawfishPath(ch_sawfish_bw_or_bedgraph) {
    ch_sawfish_bw_or_bedgraph
        .transpose()
        .map { meta, file ->
            def sample_id = file.parent.name.replaceFirst(/[^_]*_/, "")
            [meta + ['id': sample_id], file]
        }
}
