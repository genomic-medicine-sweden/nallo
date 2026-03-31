include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main'
include { WHATSHAP_HAPLOTAG } from '../../../modules/nf-core/whatshap/haplotag/main'
include { WHATSHAP_PHASE    } from '../../../modules/nf-core/whatshap/phase/main'

workflow WHATSHAP {
    take:
    ch_snv_vcf       // channel: [ val(meta), path(vcf) ]
    ch_snv_index     // channel: [ val(meta), path(tbi) ]
    ch_bam_bai       // channel: [ val(meta), path(bam), path(bai) ]
    fasta            // channel: [ val(meta), path(fasta) ]
    fai              // channel: [ val(meta), path(fai) ]
    ch_pedigree      // channel: [ val(meta), path(pedigree) ]

    main:
    // Fix metadata to group by family
    ch_bam_bai
        .map { meta, bam, bai -> [[id: meta.family_id], bam, bai] }
        .groupTuple()
        .set { ch_bam_bai_grouped }

    // Join VCFS, then join with BAMs to ensure input channel order
    // The joined VCFs and BAMs are then separated so we can pass them into WhatsHap
    ch_snv_vcf
        .map { meta, vcf -> [[id: meta.id], vcf] }
        .join(ch_snv_index, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_bam_bai_grouped, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_pedigree, failOnMismatch: true, failOnDuplicate: true)
        .multiMap { meta, snv, tbi, bam, bai, pedigree ->
            vcf: [meta, snv, tbi]
            bam: [meta, bam, bai]
            ped: [meta, pedigree]
        }
        .set { ch_whatshap_phase_in }

    fasta
        .join(fai, failOnMismatch: true, failOnDuplicate: true)
        .first()
        .set { ch_fasta_fai }

    WHATSHAP_PHASE(
        ch_whatshap_phase_in.vcf,
        ch_whatshap_phase_in.bam,
        ch_fasta_fai,
        ch_whatshap_phase_in.ped
    )

    // We cannot use the grouped BAM channel here because WhatsHap can haplotag only one BAM at a time.
    // Using combine instead of join because the VCFs are family-level, not sample-level
    // Therefore, there might be multiple BAMs per VCF and join only keeps the first match
    // (unless failOnDuplicate is true, then we get an error)

    ch_bam_bai
        .map { meta, bam, bai -> [[id: meta.family_id], meta, bam, bai] }
        .combine(WHATSHAP_PHASE.out.vcf, by: 0)
        .combine(WHATSHAP_PHASE.out.tbi, by: 0)
        .map { _family_meta, sample_meta, bam, bai, vcf, tbi -> [sample_meta, vcf, tbi, bam, bai] }
        .set { ch_whatshap_haplotag_in }

    WHATSHAP_HAPLOTAG(
        ch_whatshap_haplotag_in,
        fasta,
        fai,
        false
    )

    SAMTOOLS_INDEX(
        WHATSHAP_HAPLOTAG.out.bam
    )

    WHATSHAP_HAPLOTAG.out.bam
        .join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)
        .set { ch_bam_bai_haplotagged }

    emit:
    phased_family_snvs     = WHATSHAP_PHASE.out.vcf // channel: [ val(meta), path(vcf) ]
    phased_family_snvs_tbi = WHATSHAP_PHASE.out.tbi // channel: [ val(meta), path(tbi) ]
    haplotagged_bam_bai    = ch_bam_bai_haplotagged // channel: [ val(meta), path(bam), path(bai) ]
}
