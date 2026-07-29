include { SVDB_MERGE as SVDB_MERGE_BY_CALLER } from '../../../modules/nf-core/svdb/merge/main'
include { SVDB_MERGE as SVDB_MERGE_BY_FAMILY } from '../../../modules/nf-core/svdb/merge/main'

workflow MERGE_SV_VCFS {
    take:
    ch_reheadered_vcf   // channel: [ val(meta), path(vcf) ]
    ch_vcf              // channel: [ val(meta), path(vcf) ]
    sv_callers_to_merge // List: [ 'caller1', 'caller2', 'caller3' ]
    caller_priority     // List: [ 'caller3', 'caller1', 'caller2' ]

    main:
    // Merge the reheadered SV calls with the ones that didn't need reheadering
    ch_svdb_merge_by_caller_input = ch_reheadered_vcf
        .concat(ch_vcf)
        .map { meta, vcf -> [['id': meta.family_id, 'sv_caller': meta.sv_caller], vcf] }
        .groupTuple()

    /*
     * First merge SV calls from each caller into family VCFs
     * HiFiCNV has a different BND distance from the other callers,
     * Sawfish is not really merged (run with --no_intra), unless we are forcing joint-calling single samples and using SVDB for merging.
     * These options are set in the config
     */
    SVDB_MERGE_BY_CALLER(
        ch_svdb_merge_by_caller_input,
        [],
        true,
    )

    /*
     * Then merge the family VCFs for each caller into a single family VCF.
     * First we need to filter the SV callers to merge,
     * Then we need to group by family (meta.id), and sort the VCFs by the caller priority for SVDB merge.
     */
    ch_svdb_merge_by_family_input = SVDB_MERGE_BY_CALLER.out.vcf
        .filter { meta, _vcf ->
            sv_callers_to_merge.contains(meta.sv_caller)
        }
        .map { meta, vcf ->
            [meta - meta.subMap('sv_caller'), [meta.sv_caller, vcf]]
        }
        .groupTuple(
            sort: { a, b ->
                caller_priority.indexOf(a[0]) <=> caller_priority.indexOf(b[0])
            }
        )
        .map { meta, callers_vcfs ->
            def vcf_paths = callers_vcfs.collect { caller_vcf_pair -> caller_vcf_pair[1] }
            [meta, vcf_paths]
        }
    SVDB_MERGE_BY_FAMILY(
        ch_svdb_merge_by_family_input,
        caller_priority,
        true,
    )

    emit:
    family_caller_vcf = SVDB_MERGE_BY_CALLER.out.vcf // channel: [ val(meta), path(vcf) ]
    family_caller_tbi = SVDB_MERGE_BY_CALLER.out.tbi // channel: [ val(meta), path(tbi) ]
    family_vcf        = SVDB_MERGE_BY_FAMILY.out.vcf // channel: [ val(meta), path(vcf) ]
    family_tbi        = SVDB_MERGE_BY_FAMILY.out.tbi // channel: [ val(meta), path(tbi) ]

}
