include { BCFTOOLS_INDEX        } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_MERGE        } from '../../../modules/nf-core/bcftools/merge/main'
include { LONGPHASE_HAPLOTAG    } from '../../../modules/nf-core/longphase/haplotag/main'
include { LONGPHASE_PHASE       } from '../../../modules/nf-core/longphase/phase/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { SPLIT_MULTISAMPLE_VCF } from '../../../subworkflows/local/split_multisample_vcf/main'

workflow LONGPHASE {
    take:
    ch_snv_vcf           // channel: [ val(meta), path(vcf) ]
    ch_sv_vcf            // channel: [ val(meta), path(vcf) ] Optional
    ch_bam_bai           // channel: [ val(meta), path(bam), path(bai) ]
    ch_family_to_samples // channel: [ val(meta), val(list_of_sample_ids) ]
    fasta                // channel: [ val(meta), path(fasta) ]
    fai                  // channel: [ val(meta), path(fai) ]
    phase_with_svs       // bool: Whether to include SVs in phasing (true) or not (false)

    main:
    ch_snv_with_type = ch_snv_vcf
        .map { meta, vcf -> [meta, vcf, "snv"] }

    if (phase_with_svs) {
        ch_split_in = ch_sv_vcf
            .map { meta, vcf -> [meta, vcf, "sv"] }
            .mix(ch_snv_with_type)
    }
    else {
        ch_split_in = ch_snv_with_type
    }

    SPLIT_MULTISAMPLE_VCF(
        ch_split_in,
        ch_family_to_samples,
    )

    ch_split_vcfs = SPLIT_MULTISAMPLE_VCF.out.split_vcf
        .branch { meta, vcf, variant_type ->
            sv: variant_type == 'sv'
            [meta, vcf]
            snv: variant_type == 'snv'
            [meta, vcf]
        }

    ch_bam_vcf = ch_bam_bai
        .map { meta, bam, bai -> [[id: meta.id, family_id: meta.family_id], bam, bai] }
        .join(ch_split_vcfs.snv, failOnMismatch: true, failOnDuplicate: true)

    if (phase_with_svs) {
        ch_longphase_phase_in = ch_bam_vcf
            .join(ch_split_vcfs.sv, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, bam, bai, snvs, svs -> [meta, bam, bai, snvs, svs, []] }
    }
    else {
        ch_longphase_phase_in = ch_bam_vcf
            .map { meta, bam, bai, snvs -> [meta, bam, bai, snvs, [], []] }
    }

    LONGPHASE_PHASE(
        ch_longphase_phase_in,
        fasta,
        fai,
    )

    ch_bcftools_index_in = LONGPHASE_PHASE.out.snv_vcf
        .map { meta, vcf -> [meta + [variant_type: 'snv'], vcf] }
        .mix(LONGPHASE_PHASE.out.sv_vcf.map { meta, vcf -> [meta + [variant_type: 'sv'], vcf] })

    // Index all phased VCFs, ignoring variant types.
    BCFTOOLS_INDEX(ch_bcftools_index_in)

    ch_phased_vcf = ch_bcftools_index_in
        .join(BCFTOOLS_INDEX.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi -> [meta + [id: meta.family_id] - meta.subMap('family_id'), vcf, tbi] }
        .groupTuple()
        .map { meta, vcfs, tbis -> [meta, vcfs, tbis, []] }

    BCFTOOLS_MERGE(
        ch_phased_vcf,
        fasta.join(fai, failOnMismatch: true, failOnDuplicate: true).collect(),
    )

    ch_phased_family_vcfs = BCFTOOLS_MERGE.out.vcf
        .branch { meta, vcf ->
            snv: meta.variant_type == 'snv'
            [meta - meta.subMap('variant_type'), vcf]
            sv: meta.variant_type == 'sv'
            [meta - meta.subMap('variant_type'), vcf]
        }


    ch_phased_family_vcf_index = BCFTOOLS_MERGE.out.index
        .branch { meta, tbi ->
            snv: meta.variant_type == 'snv'
            [meta - meta.subMap('variant_type'), tbi]
            sv: meta.variant_type == 'sv'
            [meta - meta.subMap('variant_type'), tbi]
        }

    ch_phased_family_snvs = ch_phased_family_vcfs.snv
    ch_phased_family_snvs_tbi = ch_phased_family_vcf_index.snv
    ch_phased_family_svs = ch_phased_family_vcfs.sv
    ch_phased_family_svs_tbi = ch_phased_family_vcf_index.sv


    // Haplotagging
    if (phase_with_svs) {
        ch_vcfs_for_haplotag = LONGPHASE_PHASE.out.snv_vcf
            .join(LONGPHASE_PHASE.out.sv_vcf, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, snvs, svs -> [meta, snvs, svs, []] }
    }
    else {
        ch_vcfs_for_haplotag = LONGPHASE_PHASE.out.snv_vcf
            .map { meta, vcf -> [meta, vcf, [], []] }
    }

    // Making sure to keep the full meta we get in case downstream processes need it
    ch_longphase_haplotag_in = ch_bam_bai
        .map { meta, bam, bai -> [[id: meta.id, family_id: meta.family_id], meta, bam, bai] }
        .join(ch_vcfs_for_haplotag, failOnMismatch: true, failOnDuplicate: true)
        .map { _join_key, full_meta, bam, bai, vcf_snvs, vcf_svs, vcf_mods ->
            [full_meta, bam, bai, vcf_snvs, vcf_svs, vcf_mods]
        }

    LONGPHASE_HAPLOTAG(
        ch_longphase_haplotag_in,
        fasta,
        fai,
    )

    SAMTOOLS_INDEX(
        LONGPHASE_HAPLOTAG.out.bam
    )

    ch_bam_bai_haplotagged = LONGPHASE_HAPLOTAG.out.bam
        .join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)

    emit:
    phased_family_snvs = ch_phased_family_snvs         // channel: [ val(meta), path(vcf) ]
    phased_family_snvs_tbi = ch_phased_family_snvs_tbi // channel: [ val(meta), path(tbi) ]
    phased_family_svs = ch_phased_family_svs           // channel: [ val(meta), path(vcf) ]
    phased_family_svs_tbi = ch_phased_family_svs_tbi   // channel: [ val(meta), path(tbi) ]
    haplotagged_bam_bai = ch_bam_bai_haplotagged       // channel: [ val(meta), path(bam), path(bai) ]
}
