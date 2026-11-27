//
// A subworkflow to add most severe consequence and pli to a vep annotated vcf
//

include { CUSTOM_ADDMOSTSEVERECONSEQUENCE } from '../../../modules/nf-core/custom/addmostsevereconsequence/main'
include { CUSTOM_ADDMOSTSEVEREPLI         } from '../../../modules/nf-core/custom/addmostseverepli/main'
include { HTSLIB_BGZIPTABIX               } from '../../../modules/nf-core/htslib/bgziptabix/main'
workflow ANNOTATE_CSQ_PLI {
    take:
    ch_vcf                  // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_variant_consequences // channel: [mandatory] [ val(meta), path(consequences) ]

    main:
    CUSTOM_ADDMOSTSEVERECONSEQUENCE(ch_vcf, ch_variant_consequences)

    CUSTOM_ADDMOSTSEVEREPLI(CUSTOM_ADDMOSTSEVERECONSEQUENCE.out.vcf)

    HTSLIB_BGZIPTABIX(
        CUSTOM_ADDMOSTSEVEREPLI.out.vcf.map { meta, vcf -> [meta, vcf, [], []] },
        'compress',
        true,
        'vcf',
    )

    emit:
    vcf = HTSLIB_BGZIPTABIX.out.output // channel: [ val(meta), path(vcf) ]
    tbi = HTSLIB_BGZIPTABIX.out.index  // channel: [ val(meta), path(tbi) ]
}
