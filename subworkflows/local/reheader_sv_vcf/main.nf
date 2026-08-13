include { BCFTOOLS_QUERY                     } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER                  } from '../../../modules/nf-core/bcftools/reheader/main'
include { GAWK as CREATE_SAMPLES_FILE        } from '../../../modules/nf-core/gawk/main'


workflow REHEADER_SV_VCF {

    take:
        ch_vcf // channel: [ val(meta), path(vcf) ]

    main:
        ch_vcf_reheadered = channel.empty()

    // If Severus or Sniffles was used, we need to reheader the VCF
    // Since Sniffles hardcodes the sample name as SAMPLE, and Severus bases it on the file name.
    // HiFiCNV doesn't have this issue, so we filter it out here, and add it back later.

    // Branching channel to get the VCFs that need reheadering and those that don't
    ch_vcf_reheader = ch_vcf.branch { meta, _vcf ->
        def callers_needing_reheader = ['severus', 'sniffles']
        to_reheader: callers_needing_reheader.contains(meta.sv_caller)
        no_reheader: !callers_needing_reheader.contains(meta.sv_caller)
    }

    // Getting the sample name from the VCFs that need reheadering
    BCFTOOLS_QUERY(
        ch_vcf_reheader.to_reheader.map { meta, vcf -> [meta, vcf, []] },
        [],
        [],
        [],
    )
    // Then create a "vcf_sample_name meta.id" file for bcftools reheader
    CREATE_SAMPLES_FILE(BCFTOOLS_QUERY.out.output, [], false)

    ch_bcftools_reheader_input = ch_vcf_reheader.to_reheader
        .join(CREATE_SAMPLES_FILE.out.output, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, samples -> [meta, vcf, [], samples] }

    // Finally, reheader the VCF with meta.id as the sample name
    BCFTOOLS_REHEADER(
        ch_bcftools_reheader_input,
        [[], []],
    )

    // Concat the reheadered VCFs with the ones that didn't need reheadering to merge later
    ch_vcf_reheadered = BCFTOOLS_REHEADER.out.vcf
        .concat(ch_vcf_reheader.no_reheader)
        .map { meta, vcf -> [['id': meta.family_id, 'sv_caller': meta.sv_caller], vcf] }
        .groupTuple()

    emit:
        vcf_reheadered = ch_vcf_reheadered // channel: [ val(meta), path(vcf) ]

}
