//
// A subworkflow to score and rank variants.
//

include { GENMOD_ANNOTATE  } from '../../../modules/nf-core/genmod/annotate/main'
include { GENMOD_MODELS    } from '../../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE     } from '../../../modules/nf-core/genmod/score/main'
include { GENMOD_COMPOUND  } from '../../../modules/nf-core/genmod/compound/main'
include { BCFTOOLS_SORT    } from '../../../modules/nf-core/bcftools/sort/main'

workflow RANK_VARIANTS {

    take:
    ch_vcf                // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_pedfile            // channel: [mandatory] [ val(meta), path(ped) ]
    ch_reduced_penetrance // channel: [mandatory] [ val(meta), path(pentrance) ]
    ch_score_config       // channel: [mandatory] [ val(meta), path(ini) ]

    main:
    ch_versions = Channel.empty()

    GENMOD_ANNOTATE ( ch_vcf )
    ch_versions = ch_versions.mix(GENMOD_ANNOTATE.out.versions)

    GENMOD_MODELS (
        GENMOD_ANNOTATE.out.vcf,
        ch_pedfile.map { meta, ped -> ped },
        ch_reduced_penetrance.map { meta, file -> file }
    )
    ch_versions = ch_versions.mix(GENMOD_MODELS.out.versions)

    GENMOD_SCORE (
        GENMOD_MODELS.out.vcf,
        ch_pedfile.map { meta, ped -> ped },
        ch_score_config.map { meta, file -> file }
    )
    ch_versions = ch_versions.mix(GENMOD_SCORE.out.versions)

    GENMOD_COMPOUND ( GENMOD_SCORE.out.vcf )
    ch_versions = ch_versions.mix(GENMOD_COMPOUND.out.versions)

    BCFTOOLS_SORT ( GENMOD_COMPOUND.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    emit:
    vcf      = BCFTOOLS_SORT.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi      = BCFTOOLS_SORT.out.tbi // channel: [ val(meta), path(tbi) ]
    versions = ch_versions           // channel: [ path(versions.yml) ]
}
