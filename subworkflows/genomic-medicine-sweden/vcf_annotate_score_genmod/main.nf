include { BCFTOOLS_VIEW   } from '../../../modules/nf-core/bcftools/view/main'
include { GENMOD_ANNOTATE } from '../../../modules/nf-core/genmod/annotate/main'
include { GENMOD_COMPOUND } from '../../../modules/nf-core/genmod/compound/main'
include { GENMOD_MODELS   } from '../../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE    } from '../../../modules/nf-core/genmod/score/main'

workflow VCF_ANNOTATE_SCORE_GENMOD {
    take:
    ch_vcf                       // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_ped                       // channel: [optional]  [ val(meta), path(ped) ]
    ch_reduced_penetrance        // channel: [optional]  [ val(meta), path(penetrance) ]
    ch_score_config              // channel: [mandatory] [ val(meta), path(ini) ]
    val_run_score_only           // Boolean: If true, only run the scoring step

    main:
    def val_run_annotate_and_models = !val_run_score_only
    def val_run_compounds           = !val_run_score_only

    if (val_run_annotate_and_models) {
        GENMOD_ANNOTATE(
            ch_vcf
        )

        ch_genmod_models_in = GENMOD_ANNOTATE.out.vcf.join(ch_ped, failOnMismatch: true, failOnDuplicate: true)

        GENMOD_MODELS(
            ch_genmod_models_in,
            ch_reduced_penetrance.map { _meta, file -> file },
        )

        ch_vcf_for_genmod_score = GENMOD_MODELS.out.vcf
    }
    else {
        ch_vcf_for_genmod_score = ch_vcf
    }

    def ch_genmod_score_in = ch_vcf_for_genmod_score
        .join(ch_ped, failOnDuplicate: true, remainder: true)
        .map { meta, vcf, ped ->
            ped ? [meta, vcf, ped] : [meta, vcf, []]
        }
        .join(ch_score_config, failOnMismatch: true, failOnDuplicate: true)

    GENMOD_SCORE(
        ch_genmod_score_in
    )

    if (val_run_compounds) {
        GENMOD_COMPOUND(
            GENMOD_SCORE.out.vcf
        )

        ch_bcftools_view_in = GENMOD_COMPOUND.out.vcf
    }
    else {
        ch_bcftools_view_in = GENMOD_SCORE.out.vcf
    }

    // Genmod can only output an uncompressed VCF, bcftools view can be used to compress and index the output file
    BCFTOOLS_VIEW(
        ch_bcftools_view_in.map { meta, vcf -> [meta, vcf, []] },
        [],
        [],
        [],
    )

    emit:
    vcf   = BCFTOOLS_VIEW.out.vcf                            // channel: [ val(meta), path(vcf) ]
    index = BCFTOOLS_VIEW.out.tbi.mix(BCFTOOLS_VIEW.out.csi) // channel: [ val(meta), path(index) ]
}
