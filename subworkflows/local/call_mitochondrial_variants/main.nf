/*
 * Workflow to call mitochondrial variants
 */

include { DEEPVARIANT_RUNDEEPVARIANT } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { MITORSAW_HAPLOTYPE          } from '../../../modules/nf-core/mitorsaw/haplotype/main'

workflow CALL_MITOCHONDRIAL_VARIANTS {
    take:
    ch_bam_bai            // channel: [val(meta), path(bam), path(bai)]
    ch_fasta              // channel: [val(meta), path(fasta)]
    ch_fai                // channel: [val(meta), path(fai)]
    ch_par_bed            // channel: [val(meta), path(bed)]  – PAR regions (deepvariant)
    ch_mitochondrial_bed  // channel: [val(meta), path(bed)]  – mitochondrial interval (deepvariant)
    mitochondrial_caller  // string

    main:
    ch_fasta.join(ch_fai).set { ch_fasta_fai }

    if (mitochondrial_caller == "mitorsaw") {

        MITORSAW_HAPLOTYPE(
            ch_bam_bai,
            ch_fasta_fai,
            [],
            [],
        )

        ch_vcf = MITORSAW_HAPLOTYPE.out.vcf
        ch_tbi = MITORSAW_HAPLOTYPE.out.tbi

    } else if (mitochondrial_caller == "deepvariant") {

        // Broadcast the single mito BED to every sample, skip if BED is empty
        ch_bam_bai
            .combine(ch_mitochondrial_bed)
            .filter { _bam_meta, _bam, _bai, _mitochondrial_meta, bed -> bed.size() > 0 }
            .map { bam_meta, bam, bai, mitochondrial_meta, bed ->
                [bam_meta + [genome: mitochondrial_meta.genome, num_intervals: mitochondrial_meta.num_intervals], bam, bai, bed] }
            .set { deepvariant_in }

        DEEPVARIANT_RUNDEEPVARIANT(
            deepvariant_in,
            ch_fasta,
            ch_fai,
            [[], []],
            ch_par_bed,
        )

        ch_vcf = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
        ch_tbi = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_tbi

    }

    emit:
    mitochondrial_vcf = ch_vcf  // channel: [val(meta), path(vcf/gvcf)]
    mitochondrial_tbi = ch_tbi  // channel: [val(meta), path(tbi)]
}
