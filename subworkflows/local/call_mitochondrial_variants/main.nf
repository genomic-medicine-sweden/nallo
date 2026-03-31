/*
 * Workflow to call mitochondrial variants
 */

include { MITORSAW_HAPLOTYPE } from '../../../modules/nf-core/mitorsaw/haplotype/main'

workflow CALL_MITOCHONDRIAL_VARIANTS {
    take:
    ch_bam_bai
    ch_fasta_fai
    mitochondrial_caller

    main:
    if (mitochondrial_caller == "mitorsaw") {

        MITORSAW_HAPLOTYPE(
            ch_bam_bai,
            ch_fasta_fai,
        )

        ch_vcf = MITORSAW_HAPLOTYPE.out.vcf
        ch_tbi = MITORSAW_HAPLOTYPE.out.tbi

    }

    emit:
    mitochondrial_vcf = ch_vcf
    mitochondrial_tbi = ch_tbi
}
