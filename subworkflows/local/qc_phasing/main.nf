include { BCFTOOLS_CONCAT  } from '../../../modules/nf-core/bcftools/concat/main'
include { CRAMINO          } from '../../../modules/nf-core/cramino/main'
include { WHATSHAP_STATS   } from '../../../modules/nf-core/whatshap/stats/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow QC_PHASING {
    take:
    ch_phased_family_snvs     // channel: [ val(meta), path(vcf) ]
    ch_phased_family_snvs_tbi // channel: [ val(meta), path(tbi) ]
    ch_phased_family_svs      // channel: [ val(meta), path(vcf) ] Optional
    ch_phased_family_svs_tbi  // channel: [ val(meta), path(tbi) ] Optional
    ch_bam_bai_haplotagged    // channel: [ val(meta), path(bam), path(bai) ]
    ch_family_to_samples      // channel: [ val(family_id), val(list_of_sample_ids) ]
    phase_with_svs            //    bool: Whether SVs were included in phasing (true) or not (false)
    run_whatshap_stats        //    bool: Whether to run WHATSHAP_STATS (true) or not (false)

    main:
    ch_versions = channel.empty()
    ch_whatshap_stats_tsv = channel.empty()
    ch_whatshap_stats_gtf = channel.empty()

    // If we co-phased SVs, concatenate SNV and SV VCFs to get accurate stats from WhatsHap
    if (phase_with_svs) {
        ch_phased_family_snvs
            .join(ch_phased_family_snvs_tbi, failOnMismatch: true, failOnDuplicate: true)
            .mix(
                ch_phased_family_svs.join(ch_phased_family_svs_tbi, failOnMismatch: true, failOnDuplicate: true)
            )
            .groupTuple()
            .set { ch_bcftools_concat_in }

        BCFTOOLS_CONCAT(ch_bcftools_concat_in).vcf.set { ch_phased_vcf }
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
    }
    else {
        ch_phased_family_snvs.set { ch_phased_vcf }
    }

    // At this point, we have exactly one phased VCF per family.
    // Whatshap stats only works on one sample at a time, even if there are multiple samples in the VCF.
    // Therefore, we join with the known samples per family and create one item per sample,
    // duplicating the VCF and TBI paths as needed.

    ch_phased_vcf
        .join(ch_family_to_samples, failOnMismatch: true, failOnDuplicate: true)
        .transpose()
        .map { meta, vcf, sample_id ->
            [meta + [id: sample_id, family_id: meta.id], vcf]
        }
        .set { ch_phased_vcf }

    if(run_whatshap_stats) {

        WHATSHAP_STATS(
            ch_phased_vcf,
            true,
            true,
            false,
        )

        ch_whatshap_stats_tsv = WHATSHAP_STATS.out.tsv
        ch_whatshap_stats_gtf = WHATSHAP_STATS.out.gtf

    }

    TABIX_BGZIPTABIX(ch_whatshap_stats_gtf)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    TABIX_BGZIPTABIX.out.gz_tbi
        .multiMap { meta, gtf_gz, gtf_tbi ->
            gz: [meta, gtf_gz]
            tbi: [meta, gtf_tbi]
        }
        .set { ch_phasing_gtf }

    CRAMINO(ch_bam_bai_haplotagged)
    ch_versions = ch_versions.mix(CRAMINO.out.versions)

    emit:
    phasing_stats        = ch_phasing_stats       // channel: [ val(meta), path(stats) ]
    phasing_blocks       = ch_phasing_gtf.gz      // channel: [ val(meta), path(gtf) ]
    phasing_blocks_index = ch_phasing_gtf.tbi     // channel: [ val(meta), path(gtf_index) ]
    haplotagging_stats   = CRAMINO.out.stats      // channel: [ val(meta), path(stats) ]
    versions             = ch_versions            // channel: [ path(versions.yml) ]
}
