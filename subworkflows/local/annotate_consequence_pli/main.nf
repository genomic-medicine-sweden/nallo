//
// A subworkflow to add most severe consequence and pli to a vep annotated vcf
//

include { CUSTOM_ADDMOSTSEVERECONSEQUENCE } from '../../../modules/nf-core/custom/addmostsevereconsequence/main'
include { CUSTOM_ADDMOSTSEVEREPLI         } from '../../../modules/nf-core/custom/addmostseverepli/main'
include { TABIX_TABIX                     } from '../../../modules/nf-core/tabix/tabix/main'
workflow ANNOTATE_CSQ_PLI {
    take:
    ch_vcf                  // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_variant_consequences // channel: [mandatory] [ val(meta), path(consequences) ]

    main:
    CUSTOM_ADDMOSTSEVERECONSEQUENCE (ch_vcf, ch_variant_consequences)

    CUSTOM_ADDMOSTSEVEREPLI (CUSTOM_ADDMOSTSEVERECONSEQUENCE.out.vcf)

    TABIX_TABIX (CUSTOM_ADDMOSTSEVEREPLI.out.vcf)

    emit:
    vcf = CUSTOM_ADDMOSTSEVEREPLI.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi = TABIX_TABIX.out.index           // channel: [ val(meta), path(tbi) ]
}
