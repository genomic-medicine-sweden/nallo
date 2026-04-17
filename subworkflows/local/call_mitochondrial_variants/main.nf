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

        // Broadcast the single mito BED to every sample
        ch_bam_bai
            .dump(tag: "bam_bai in mitochondrial workflow")
            .combine(ch_mitochondrial_bed)
            .dump(tag: "bam_bai input to call_snvs mitochondrial workflow")
            .map { bam_meta, bam, bai, mito_meta, bed ->
                [bam_meta + [genome: mito_meta.genome, num_intervals: 1], bam, bai, bed, 1] }
            .set { call_snvs_input }
        call_snvs_input.dump(tag: "input to call_snvs mitochondrial workflow")

        DEEPVARIANT_RUNDEEPVARIANT(
            call_snvs_input,
            ch_fasta,
            ch_fai,
            [[], []],
            ch_par_bed,
        )

        ch_vcf = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
        ch_tbi = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_tbi

    }

    emit:
    mitochondrial_vcf = ch_vcf  // channel: [val(meta), path(vcf/gvcf)]
    mitochondrial_tbi = ch_tbi  // channel: [val(meta), path(tbi)]
}
