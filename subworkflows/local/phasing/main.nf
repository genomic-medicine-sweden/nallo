include { HIPHASE          } from '../../../subworkflows/local/hiphase/main'
include { LONGPHASE        } from '../../../subworkflows/local/longphase/main'
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'
include { WHATSHAP         } from '../../../subworkflows/local/whatshap/main'
include { QC_PHASING       } from '../../../subworkflows/local/qc_phasing/main'


workflow PHASING {
    take:
    ch_snv_vcf // channel: [ val(meta), path(vcf) ]
    ch_snv_vcf_index // channel: [ val(meta), path(tbi) ]
    ch_sv_vcf // channel: [ val(meta), path(vcf) ] Optional
    ch_sv_vcf_index // channel: [ val(meta), path(tbi) ] Optional
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    ch_family_to_samples // channel: [ val(meta), val(set_of_sample_ids) ]
    fasta // channel: [ val(meta), path(fasta) ]
    fai // channel: [ val(meta), path(fai) ]
    phaser // string:  Phasing tool to use
    phase_with_svs // bool:    Whether to include SVs in phasing (true) or not (false)
    cram_output // bool:    Publish alignments as CRAM (true) or BAM (false)

    main:
    // Phase variants and haplotag reads with Longphase
    if (phaser.equals("longphase")) {

        LONGPHASE(
            ch_snv_vcf,
            ch_sv_vcf,
            ch_bam_bai,
            ch_family_to_samples,
            fasta,
            fai,
            phase_with_svs,
        )

        ch_phased_family_snvs = LONGPHASE.out.phased_family_snvs
        ch_phased_family_snvs_tbi = LONGPHASE.out.phased_family_snvs_tbi
        ch_phased_family_svs = phase_with_svs ? LONGPHASE.out.phased_family_svs : ch_sv_vcf
        ch_phased_family_svs_tbi = phase_with_svs ? LONGPHASE.out.phased_family_svs_tbi : ch_sv_vcf_index
        ch_bam_bai_haplotagged = LONGPHASE.out.haplotagged_bam_bai
    // Phase variants and haplotag reads with WhatsHap
    } else if (phaser.equals("whatshap")) {

        WHATSHAP(
            ch_snv_vcf,
            ch_snv_vcf_index,
            ch_bam_bai,
            fasta,
            fai,
        )

        ch_phased_family_snvs = WHATSHAP.out.phased_family_snvs
        ch_phased_family_snvs_tbi = WHATSHAP.out.phased_family_snvs_tbi
        ch_phased_family_svs = ch_sv_vcf
        ch_phased_family_svs_tbi = ch_sv_vcf_index
        ch_bam_bai_haplotagged = WHATSHAP.out.haplotagged_bam_bai
    // Phase variants and haplotag reads with HiPhase
    } else if (phaser.equals("hiphase")) {

        HIPHASE(
            ch_snv_vcf,
            ch_snv_vcf_index,
            ch_sv_vcf,
            ch_sv_vcf_index,
            ch_bam_bai,
            ch_family_to_samples,
            fasta,
            fai,
            phase_with_svs,
        )

        ch_phased_family_snvs = HIPHASE.out.phased_snvs
        ch_phased_family_snvs_tbi = HIPHASE.out.phased_snvs_tbi
        ch_phased_family_svs = HIPHASE.out.phased_svs
        ch_phased_family_svs_tbi = HIPHASE.out.phased_svs_tbi
        ch_bam_bai_haplotagged = HIPHASE.out.haplotagged_bam_bai
    }

    QC_PHASING(
        ch_phased_family_snvs,
        ch_phased_family_snvs_tbi,
        ch_phased_family_svs,
        ch_phased_family_svs_tbi,
        ch_bam_bai_haplotagged,
        ch_family_to_samples,
        phase_with_svs && !phaser.equals("whatshap"),
    )

    if (cram_output) {
        SAMTOOLS_CONVERT(
            ch_bam_bai_haplotagged,
            fasta.join(fai).collect(),
        )

        ch_haplotagged_cram_crai = SAMTOOLS_CONVERT.out.cram.join(SAMTOOLS_CONVERT.out.crai)
    }

    emit:
    phased_family_snvs     = ch_phased_family_snvs // channel: [ val(meta), path(vcf) ]
    phased_family_snvs_tbi = ch_phased_family_snvs_tbi // channel: [ val(meta), path(tbi) ]
    phased_family_svs      = ch_phased_family_svs // channel: [ val(meta), path(vcf) ]
    phased_family_svs_tbi  = ch_phased_family_svs_tbi // channel: [ val(meta), path(tbi) ]
    haplotagged_bam_bai    = ch_bam_bai_haplotagged // channel: [ val(meta), path(bam), path(bai) ]
    haplotagged_cram_crai  = cram_output ? ch_haplotagged_cram_crai : channel.empty() // channel: [ val(meta), path(cram), path(crai) ]
    stats                  = QC_PHASING.out.phasing_stats // channel: [ val(meta), path("*.stats.tsv") ]
    blocks                 = QC_PHASING.out.phasing_blocks // channel: [ val(meta), path("*.blocks.gtf.gz") ]
    blocks_index           = QC_PHASING.out.phasing_blocks_index // channel: [ val(meta), path("*.blocks.gtf.gz.tbi") ]
    haplotagging_stats     = QC_PHASING.out.haplotagging_stats   // channel: [ val(meta), path("*.txt") ]
    haplotagging_arrow     = QC_PHASING.out.haplotagging_arrow   // channel: [ val(meta), path("*.arrow") ]
}
