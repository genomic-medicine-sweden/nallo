include { HIPHASE as RUN_HIPHASE } from '../../../modules/nf-core/hiphase/main'

workflow HIPHASE {
    take:
    ch_snv_vcf // channel: [ val(meta), path(vcf) ]
    ch_snv_vcf_index // channel: [ val(meta), path(tbi) ]
    ch_sv_vcf // channel: [ val(meta), path(vcf) ] Optional
    ch_sv_vcf_index // channel: [ val(meta), path(tbi) ] Optional
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    ch_family_to_samples // channel: [ val(meta), val( [ sample_ids ] ) ]
    fasta // channel: [ val(meta), path(fasta) ]
    fai // channel: [ val(meta), path(fai) ]
    phase_with_svs // bool: Whether to include SVs in phasing (true) or not (false)

    main:
    // Prepare SNV VCF with index
    ch_snv_vcf_tbi = ch_snv_vcf.join(ch_snv_vcf_index, failOnMismatch: true, failOnDuplicate: true)

    // Group BAM files by family and join with SNV VCF
    ch_hiphase_bam_snv = ch_bam_bai
        .map { meta, bam, bai -> [[id: meta.family_id], bam, bai] }
        .groupTuple()
        .join(ch_snv_vcf_tbi, failOnMismatch: true, failOnDuplicate: true)

    // Prepare input based on whether SVs are included
    if (phase_with_svs) {
        ch_sv_vcf_tbi = ch_sv_vcf.join(ch_sv_vcf_index, failOnMismatch: true, failOnDuplicate: true)

        ch_bam_vcf = ch_hiphase_bam_snv.join(ch_sv_vcf_tbi, failOnMismatch: true, failOnDuplicate: true)
    }
    else {
        ch_bam_vcf = ch_hiphase_bam_snv.map { meta, bams, bais, snv_vcf, snv_tbi -> [meta, bams, bais, snv_vcf, snv_tbi, [], []] }
    }

    // Adding sample IDs to input tuples
    ch_hiphase_in = ch_bam_vcf.join(ch_family_to_samples, failOnMismatch: true, failOnDuplicate: true)
    // Run HiPhase
    RUN_HIPHASE(
        ch_hiphase_in,
        fasta.join(fai, failOnMismatch: true, failOnDuplicate: true).collect(),
        true,
        false,
        false,
        false,
        false,
        'tsv',
    )

    // Prepare haplotagged BAM output by matching with original metadata
    // The HiPhase output channels only contain family IDs
    // We need to match the haplotagged BAMs back to the original sample metadata
    // as to not lose information downstream processes might depend on.
    ch_haplotagged_bam_bai = RUN_HIPHASE.out.bams
        .join(RUN_HIPHASE.out.bams_indexes, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, bams, bais ->
            // Keep BAM and BAI pairing while enforcing deterministic order.
            def bam_bai_pairs = [bams, bais]
                .transpose()
                .sort { left, right ->
                    left[0].getName() <=> right[0].getName()
                }
            [meta, bam_bai_pairs.collect { pair -> pair[0] }, bam_bai_pairs.collect { pair -> pair[1] }]
        }
        .transpose()
        .combine(ch_bam_bai)
        .filter { _meta_phased, bam_phased, _bai_phased, meta_orig, _bam_orig, _bai_orig ->
            bam_phased.simpleName.startsWith(meta_orig.id)
        }
        .map { _meta_phased, bam_phased, bai_phased, meta_orig, _bam_orig, _bai_orig ->
            [meta_orig, bam_phased, bai_phased]
        }

    emit:
    phased_snvs         = RUN_HIPHASE.out.vcfs // channel: [ val(meta), path(vcf) ]
    phased_snvs_tbi     = RUN_HIPHASE.out.vcfs_indexes // channel: [ val(meta), path(tbi) ]
    phased_svs          = phase_with_svs ? RUN_HIPHASE.out.sv_vcfs : ch_sv_vcf // channel: [ val(meta), path(vcf) ]
    phased_svs_tbi      = phase_with_svs ? RUN_HIPHASE.out.sv_vcfs_indexes : ch_sv_vcf_index // channel: [ val(meta), path(tbi) ]
    haplotagged_bam_bai = ch_haplotagged_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
}
