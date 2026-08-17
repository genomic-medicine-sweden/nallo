include { BCFTOOLS_QUERY              } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER           } from '../../../modules/nf-core/bcftools/reheader/main'
include { GAWK as CREATE_SAMPLES_FILE } from '../../../modules/nf-core/gawk/main'


workflow REHEADER_SV_VCF {
    take:
    ch_vcf // channel: [ val(meta), path(vcf) ]

    main:
    // Sniffles hardcodes the sample name as SAMPLE, and Severus bases it on the file name,
    // so those need reheadering. HiFiCNV and sawfish don't have this issue.
    def caller_needs_reheader = [
        'severus' : true,
        'sniffles': true,
        'hificnv' : false,
        'sawfish' : false,
    ]

    // Branching channel to get the VCFs that need reheadering and those that don't. If the sv caller is not in the map keyset, throw an error.
    ch_vcf_reheader = ch_vcf.branch { meta, _vcf ->
        if (!(meta.sv_caller in caller_needs_reheader.keySet())) {
            error(
                "Unknown sv_caller '${meta.sv_caller}' in REHEADER_SV_VCF. " +
                "Allowed values: ${caller_needs_reheader.keySet().sort()}."
            )
        }

        def needs_reheader = caller_needs_reheader[meta.sv_caller]
        reheader: needs_reheader
        no_reheader: !needs_reheader
        }

    // Getting the sample name from the VCFs that need reheadering
    BCFTOOLS_QUERY(
        ch_vcf_reheader.reheader.map { meta, vcf -> [meta, vcf, []] },
        [],
        [],
        [],
    )

    // Then create a "vcf_sample_name meta.id" file for bcftools reheader
    CREATE_SAMPLES_FILE(BCFTOOLS_QUERY.out.output, [], false)

    ch_bcftools_reheader_input = ch_vcf_reheader.reheader
        .join(CREATE_SAMPLES_FILE.out.output, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, samples -> [meta, vcf, [], samples] }

    // Finally, reheader the VCF with meta.id as the sample name
    BCFTOOLS_REHEADER(
        ch_bcftools_reheader_input,
        [[], []],
    )

    // Concat the reheadered VCFs and indices with the ones that didn't need reheadering to merge later
    ch_vcf_reheadered = BCFTOOLS_REHEADER.out.vcf
        .concat(ch_vcf_reheader.no_reheader)

    ch_tbi_reheadered = BCFTOOLS_REHEADER.out.index
        .concat(ch_vcf_reheader.no_reheader.map { meta, vcf -> [meta, vcf + '.tbi'] })

    emit:
    vcf = ch_vcf_reheadered // channel: [ val(meta), path(vcf) ]
    tbi = ch_tbi_reheadered // channel: [ val(meta), path(tbi) ]
}
