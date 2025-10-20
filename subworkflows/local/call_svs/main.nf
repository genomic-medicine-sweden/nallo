include { ADD_FOUND_IN_TAG                          } from '../../../modules/local/add_found_in_tag/main'
include { CLEAN_SNIFFLES                            } from '../../../modules/local/clean_sniffles/main'
include { SVDB_MERGE as SVDB_MERGE_BY_CALLER        } from '../../../modules/nf-core/svdb/merge/main'
include { SVDB_MERGE as SVDB_MERGE_BY_FAMILY        } from '../../../modules/nf-core/svdb/merge/main'
include { BCFTOOLS_VIEW                             } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_QUERY                            } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER                         } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_SORT                             } from '../../../modules/nf-core/bcftools/sort/main'
include { CREATE_SAMPLES_FILE                       } from '../../../modules/local/create_samples_file/main'
include { HIFICNV                                   } from '../../../modules/local/pacbio/hificnv'
include { SAWFISH_DISCOVER                          } from '../../../modules/nf-core/sawfish/discover/main'
include { SAWFISH_JOINTCALL                         } from '../../../modules/nf-core/sawfish/jointcall/main'
include { SEVERUS                                   } from '../../../modules/nf-core/severus/main'
include { SNIFFLES                                  } from '../../../modules/nf-core/sniffles/main'
include { TABIX_TABIX as TABIX_HIFICNV              } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as TABIX_SEVERUS         } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow CALL_SVS {

    take:
    ch_bam_bai                              // channel: [ val(meta), path(bam), path(bai) ]
    ch_tandem_repeats                       // channel: [ val(meta), path(bed) ]
    ch_snvs                                 // channel: [ val(meta), path(vcf) ]
    ch_fasta                                // channel: [ val(meta), path(fasta) ]
    ch_expected_xy_bed                      // channel: [ val(meta), path(bed) ]
    ch_expected_xx_bed                      // channel: [ val(meta), path(bed) ]
    ch_exclude_bed                          // channel: [ val(meta), path(bed) ]
    sv_callers_to_run                       //    List: [ 'caller1', 'caller2', 'caller3' ]
    sv_callers_to_merge                     //    List: [ 'caller1', 'caller2', 'caller3' ]
    caller_priority                         //    List: [ 'caller3', 'caller1', 'caller2' ]
    ch_sv_call_regions                      // channel: [ val(meta), path(bed) ]
    filter_calls_on_regions                 //    bool: Should we filter SV calls to the regions provided in ch_sv_call_regions?
    force_sawfish_joint_call_single_samples //    bool: Force joint-calling with Sawfish even for single samples

    main:
    ch_versions = Channel.empty()
    ch_sv_calls = Channel.empty()

    //
    // Call SVs with Severus
    //
    if(sv_callers_to_run.contains('severus')) {

        SEVERUS (
            ch_bam_bai.map { meta, bam, bai -> [ meta, bam, bai, [], [], [] ] },
            ch_tandem_repeats
        )
        ch_versions = ch_versions.mix(SEVERUS.out.versions)

        TABIX_SEVERUS (
            SEVERUS.out.all_vcf
        )
        ch_versions = ch_versions.mix(TABIX_SEVERUS.out.versions)

        ch_sv_calls = ch_sv_calls.mix(
            addCallerToMeta(
                TABIX_SEVERUS.out.gz_tbi,
                'severus'
            )
        )
    }

    //
    // Call SVs with Sniffles
    //
    if(sv_callers_to_run.contains('sniffles')) {

        SNIFFLES (
            ch_bam_bai
        )
        ch_versions = ch_versions.mix(SNIFFLES.out.versions)

        CLEAN_SNIFFLES (
            SNIFFLES.out.vcf
        )
        ch_versions = ch_versions.mix(CLEAN_SNIFFLES.out.versions)

        BCFTOOLS_SORT (
            CLEAN_SNIFFLES.out.vcf
        )
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

        ch_sv_calls = ch_sv_calls.mix(
            addCallerToMeta(
                BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_SORT.out.tbi, failOnMismatch:true, failOnDuplicate:true),
                'sniffles'
            )
        )
    }

    //
    // Call CNVs with HiFiCNV
    //
    if(sv_callers_to_run.contains('hificnv')) {

        ch_bam_bai
            .join(ch_snvs, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, bam, bai, vcf -> [ meta, bam, bai, vcf, meta.sex ] }
            .set { ch_hificnv_input }

        HIFICNV (
            ch_hificnv_input,
            ch_fasta,
            ch_expected_xy_bed,
            ch_expected_xx_bed,
            ch_exclude_bed
        )
        ch_versions = ch_versions.mix(HIFICNV.out.versions)

        TABIX_HIFICNV (
            HIFICNV.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_HIFICNV.out.versions)

        ch_sv_calls = ch_sv_calls.mix(
            addCallerToMeta(
                HIFICNV.out.vcf.join(TABIX_HIFICNV.out.tbi, failOnMismatch:true, failOnDuplicate:true),
                'hificnv'
            )
        )
    }

    //
    // Call SVs with Sawfish
    //
    if(sv_callers_to_run.contains('sawfish')) {

        ch_bam_bai
            .join(ch_snvs, failOnMismatch:true, failOnDuplicate:true)
            .combine(ch_expected_xx_bed)
            .combine(ch_expected_xy_bed)
            .multiMap { meta, bam, bai, vcf, xx_meta, xx_bed, xy_meta, xy_bed ->
                bam_bai: [ meta, bam, bai ]
                vcf: [ meta, vcf ]
                expected_copynumber_bed: meta.sex == 1
                    ? [ xy_meta, xy_bed ]
                    : [ xx_meta, xx_bed ]
            }
            .set { ch_sawfish_discover_input }

        SAWFISH_DISCOVER (
            ch_sawfish_discover_input.bam_bai,
            ch_fasta,
            ch_sawfish_discover_input.expected_copynumber_bed,
            ch_sawfish_discover_input.vcf,
            ch_exclude_bed
        )
        ch_versions = ch_versions.mix(SAWFISH_DISCOVER.out.versions)

        // Sawfish needs joint-calling to actually produce SV calls. Without it, there are no sample names
        // in the VCFs, and they can't be post-processed with bcftools. Therefore, we do joint-calling step
        // here directly, and skip doing it later with SVDB merging.
        SAWFISH_DISCOVER.out.discover_dir
            .join(ch_sawfish_discover_input.bam_bai, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, discover_dir, bam, bai ->
                def new_meta = force_sawfish_joint_call_single_samples
                    ? meta
                    : [ id: meta.family_id, family_id: meta.family_id ]

                [ new_meta, discover_dir, bam, bai ]
            }
            .groupTuple()
            .multiMap { meta, discover_dirs, bams, bais ->
                dir: [ meta, discover_dirs ]
                bam_bai: [ meta, bams, bais ]
            }
            .set { ch_sawfish_jointcall_input }

        SAWFISH_JOINTCALL (
            ch_sawfish_jointcall_input.dir,
            ch_fasta,
            ch_sawfish_jointcall_input.bam_bai,
            [[],[]]

        )
        ch_versions = ch_versions.mix(SAWFISH_JOINTCALL.out.versions)

        ch_sv_calls = ch_sv_calls.mix(
            addCallerToMeta(
                SAWFISH_JOINTCALL.out.vcf
                    .join(SAWFISH_JOINTCALL.out.tbi, failOnMismatch:true, failOnDuplicate:true),
                'sawfish'
            )
        )
    }

    //
    // Post-process SV calls
    //
    if ( filter_calls_on_regions ) {

        BCFTOOLS_VIEW (
            ch_sv_calls,
            ch_sv_call_regions.map { _meta, bed -> bed },
            [],
            []
        )
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

        ch_sv_calls_filtered = BCFTOOLS_VIEW.out.vcf
            .join(BCFTOOLS_VIEW.out.tbi, failOnMismatch:true, failOnDuplicate:true)

    } else {
        ch_sv_calls_filtered = ch_sv_calls
    }

    ch_sv_calls_filtered
        .multiMap { meta, vcf, tbi ->
            vcf: [ meta, vcf, tbi ]
            sv_caller: meta.sv_caller
        }
        .set { ch_add_found_in_tag_input }

    // Annotate with FOUND_IN tag
    ADD_FOUND_IN_TAG (
        ch_add_found_in_tag_input.vcf,
        ch_add_found_in_tag_input.sv_caller
    )
    ch_versions = ch_versions.mix(ADD_FOUND_IN_TAG.out.versions)

    // If Severus or Sniffles was used, we need to reheader the VCF
    // Since Sniffles hardcodes the sample name as SAMPLE, and Severus bases it on the file name.
    // HiFiCNV doesn't have this issue, so we filter it out here, and add it back later.

    // Starting with getting the sample name from the VCF
    ADD_FOUND_IN_TAG.out.vcf
        .join(ADD_FOUND_IN_TAG.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .branch { meta, _vcf, _tbi ->
            def callers_needing_reheader = [ 'severus', 'sniffles' ]
            to_reheader: callers_needing_reheader.contains(meta.sv_caller)
            no_reheader: !callers_needing_reheader.contains(meta.sv_caller)
        }
        .set { ch_found_in_tagged_vcf }

    BCFTOOLS_QUERY (
        ch_found_in_tagged_vcf.to_reheader,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions)

    // Then create a "vcf_sample_name meta.id" file for bcftools reheader
    CREATE_SAMPLES_FILE ( BCFTOOLS_QUERY.out.output )
    ch_versions = ch_versions.mix(CREATE_SAMPLES_FILE.out.versions)

    ch_found_in_tagged_vcf.to_reheader
        .join( CREATE_SAMPLES_FILE.out.samples, failOnMismatch:true, failOnDuplicate:true )
        .map { meta, vcf, _index, samples -> [ meta, vcf, [], samples ] }
        .set { ch_bcftools_reheader_input }

    // Finally, reheader the VCF with meta.id as the sample name
    BCFTOOLS_REHEADER (
        ch_bcftools_reheader_input,
        [[],[]]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    // Merge the reheadered SV calls with the ones that didn't need reheadering
    BCFTOOLS_REHEADER.out.vcf
        .join ( BCFTOOLS_REHEADER.out.index, failOnMismatch:true, failOnDuplicate:true )
        .concat ( ch_found_in_tagged_vcf.no_reheader )
        .map { meta, vcf, _tbi -> [ [ 'id': meta.family_id, 'sv_caller': meta.sv_caller ], vcf ] }
        .set { reheadered_vcfs }

    reheadered_vcfs
        .groupTuple()
        .set { ch_svdb_merge_by_caller_input }

    // First merge SV calls from each caller into family VCFs
    // HiFiCNV has a different BND distance from the other callers,
    // Sawfish is not really merged (run with --no_intra), unless we are force joint-calling single samples and using SVDB for merging.
    // These options are set in the config-
    SVDB_MERGE_BY_CALLER (
        ch_svdb_merge_by_caller_input,
        [],
        true
    )
    ch_versions = ch_versions.mix(SVDB_MERGE_BY_CALLER.out.versions)

    // Then merge the family VCFs for each caller into a single family VCF.
    // First we need to filter the SV callers to merge,
    // Then we need to group by family (meta.id), and sort the VCFs by the caller priority for SVDB merge.
    SVDB_MERGE_BY_CALLER.out.vcf
        .filter { meta, _vcf ->
            sv_callers_to_merge.contains(meta.sv_caller)
        }
        .map { meta, vcf ->
            [ meta - meta.subMap('sv_caller'), [ meta.sv_caller, vcf ] ]
        }
        .groupTuple(
            sort: { a, b ->
                caller_priority.indexOf(a[0]) <=> caller_priority.indexOf(b[0]) }
        )
        .map { meta, caller_vcf ->
            def vcf_paths = caller_vcf.collect { it[1] }
            [ meta, vcf_paths ]
        }
        .set { ch_svdb_merge_by_family_input }

    SVDB_MERGE_BY_FAMILY (
        ch_svdb_merge_by_family_input,
        caller_priority,
        true
    )
    ch_versions = ch_versions.mix(SVDB_MERGE_BY_FAMILY.out.versions)

    emit:
    family_caller_vcf = SVDB_MERGE_BY_CALLER.out.vcf // channel: [ val(meta), path(vcf) ]
    family_caller_tbi = SVDB_MERGE_BY_CALLER.out.tbi // channel: [ val(meta), path(tbi) ]
    family_vcf        = SVDB_MERGE_BY_FAMILY.out.vcf // channel: [ val(meta), path(vcf) ]
    family_tbi        = SVDB_MERGE_BY_FAMILY.out.tbi // channel: [ val(meta), path(tbi) ]
    versions          = ch_versions                  // channel: [ path(versions.yml) ]

}

def addCallerToMeta(ch_caller_calls, sv_caller) {
    ch_caller_calls.map { meta, vcf, tbi ->
        [ meta + [ sv_caller: sv_caller ], vcf, tbi ]
    }
}
