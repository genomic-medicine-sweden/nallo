//
// A subworkflow to score and rank variants.
//

include { GENMOD_ANNOTATE  } from '../../../modules/nf-core/genmod/annotate/main'
include { GENMOD_MODELS    } from '../../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE     } from '../../../modules/nf-core/genmod/score/main'
include { GENMOD_COMPOUND  } from '../../../modules/nf-core/genmod/compound/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'

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

    TABIX_BGZIPTABIX ( GENMOD_COMPOUND.out.vcf )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
    vcf      = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, vcf, tbi -> [ meta, vcf ] } // channel: [ val(meta), path(vcf) ]
    tbi      = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, vcf, tbi -> [ meta, tbi ] } // channel: [ val(meta), path(tbi) ]
    versions = ch_versions                                                         // channel: [ path(versions.yml) ]
}
