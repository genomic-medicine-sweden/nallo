include { BCFTOOLS_CONCAT  } from '../../../modules/nf-core/bcftools/concat/main'
include { CRAMINO          } from '../../../modules/nf-core/cramino/main'
include { WHATSHAP_STATS   } from '../../../modules/local/whatshap/stats/main'

workflow QC_PHASING {
    take:
    ch_phased_family_snvs     // channel: [ val(meta), path(vcf) ]
    ch_phased_family_snvs_tbi // channel: [ val(meta), path(tbi) ]
    ch_phased_family_svs      // channel: [ val(meta), path(vcf) ] Optional
    ch_phased_family_svs_tbi  // channel: [ val(meta), path(tbi) ] Optional
    ch_bam_bai_haplotagged    // channel: [ val(meta), path(bam), path(bai) ]
    ch_family_to_samples      // channel: [ val(family_id), val(list_of_sample_ids) ]
    phase_with_svs            // bool: Whether SVs were included in phasing (true) or not (false)

    main:
    ch_versions = channel.empty()

    // If we co-phased SVs, concatenate SNV and SV VCFs to get accurate stats from WhatsHap
    if (phase_with_svs) {
        ch_phased_family_snvs
            .join(ch_phased_family_snvs_tbi, failOnMismatch: true, failOnDuplicate: true)
            .mix(ch_phased_family_svs
                .join(ch_phased_family_svs_tbi, failOnMismatch: true, failOnDuplicate: true)
            )
            .groupTuple()
            .set { ch_bcftools_concat_in }

        BCFTOOLS_CONCAT( ch_bcftools_concat_in )
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

        BCFTOOLS_CONCAT.out.vcf
            .join( BCFTOOLS_CONCAT.out.tbi )
            .set { ch_phased_vcf_index }
    } else {
        ch_phased_family_snvs
            .join(ch_phased_family_snvs_tbi, failOnMismatch: true, failOnDuplicate: true)
            .set { ch_phased_vcf_index }
    }

    // At this point, we have exactly one phased VCF per family.
    // Whatshap stats only works on one sample at a time, even if there are multiple samples in the VCF.
    // Therefore, we join with the known samples per family and create one item per sample,
    // duplicating the VCF and TBI paths as needed.

    ch_phased_vcf_index
        .join( ch_family_to_samples, failOnMismatch: true, failOnDuplicate: true )
        .transpose() // go from set of samples to one sample per item, duplicating the rest
        .map { meta, vcf, tbi, sample_id ->
            [ meta + [ id: sample_id, family_id: meta.id ], vcf, tbi ]
        }
        .set { ch_phased_vcf_index }


    WHATSHAP_STATS ( ch_phased_vcf_index )
    ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions)

    CRAMINO ( ch_bam_bai_haplotagged )
    ch_versions = ch_versions.mix(CRAMINO.out.versions)

    emit:
    phasing_stats          = WHATSHAP_STATS.out.stats
    phasing_blocks         = WHATSHAP_STATS.out.blocks
    phasing_blocks_index   = WHATSHAP_STATS.out.blocks_index
    haplotagging_stats     = CRAMINO.out.stats

}
