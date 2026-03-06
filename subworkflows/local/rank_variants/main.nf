//
// A subworkflow to score and rank variants.
//

include { GENMOD_ANNOTATE } from '../../../modules/nf-core/genmod/annotate/main'
include { GENMOD_MODELS   } from '../../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE    } from '../../../modules/nf-core/genmod/score/main'
include { GENMOD_COMPOUND } from '../../../modules/nf-core/genmod/compound/main'
include { BCFTOOLS_SORT   } from '../../../modules/nf-core/bcftools/sort/main'

workflow RANK_VARIANTS {
    take:
    ch_vcf                       // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_ped                       // channel: [mandatory] [ val(meta), path(ped) ]
    ch_genmod_reduced_penetrance // channel: [mandatory] [ val(meta), path(penetrance) ]
    ch_score_config              // channel: [mandatory] [ val(meta), path(ini) ]

    main:
    GENMOD_ANNOTATE(
        ch_vcf
    )

    GENMOD_ANNOTATE.out.vcf
        .join(ch_ped, failOnMismatch: true, failOnDuplicate: true)
        .set { genmod_models_in }

    GENMOD_MODELS(
        genmod_models_in,
        ch_genmod_reduced_penetrance.map { _meta, file -> file },
    )

    GENMOD_MODELS.out.vcf
        .join(ch_ped, failOnMismatch: true, failOnDuplicate: true)
        .set { genmod_score_in }

    GENMOD_SCORE(
        genmod_score_in,
        ch_score_config.map { _meta, file -> file },
    )

    GENMOD_COMPOUND(
        GENMOD_SCORE.out.vcf
    )

    BCFTOOLS_SORT(
        GENMOD_COMPOUND.out.vcf
    )

    emit:
    vcf      = BCFTOOLS_SORT.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi      = BCFTOOLS_SORT.out.tbi // channel: [ val(meta), path(tbi) ]
}
