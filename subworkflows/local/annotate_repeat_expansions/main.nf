include { STRANGER     } from '../../../modules/nf-core/stranger/'
include { STRDROP_CALL } from '../../../modules/nf-core/strdrop/call/main'
workflow ANNOTATE_REPEAT_EXPANSIONS {
    take:
    vcf                       // channel: [ val(meta), path(vcf) ]
    strdrop_training_set_json // channel: [ val(meta), path(json) ]
    strdrop_training_set_vcfs // channel: [ val(meta), [ path(vcf) ] ]
    stranger_variant_catalog  // channel: [ val(meta), path(json) ]

    main:
    STRDROP_CALL(
        vcf,
        strdrop_training_set_json,
        strdrop_training_set_vcfs,
    )

    STRANGER(
        STRDROP_CALL.out.vcf,
        stranger_variant_catalog,
    )

    emit:
    vcf = STRANGER.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi = STRANGER.out.tbi // channel: [ val(meta), path(tbi) ]
}
