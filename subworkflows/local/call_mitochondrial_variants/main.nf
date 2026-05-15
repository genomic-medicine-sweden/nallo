/*
 * Workflow to call mitochondrial variants
 */

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_MITO } from '../../../modules/nf-core/bcftools/view/main'
include { DEEPVARIANT_RUNDEEPVARIANT          } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { MITORSAW_HAPLOTYPE                  } from '../../../modules/nf-core/mitorsaw/haplotype/main'

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

        // Add the mitochondrial BED to every sample, skip if BED is empty
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

    // Split VCF into SNVs/small indels and SVs for callers that produce both. The logic is in the config.
    // Deepvariant is SNV-only, so no split is needed.
    if (mitochondrial_caller != "deepvariant") {

        ch_vcf
            .map { meta, vcf -> [meta + [variant_type: "snv"], vcf, []] }
            .mix(ch_vcf.map { meta, vcf -> [meta + [variant_type: "sv"], vcf, []] })
            .set { ch_mito_split_input }

        BCFTOOLS_VIEW_MITO(ch_mito_split_input, [], [], [])

        BCFTOOLS_VIEW_MITO.out.vcf
            .branch { meta, _vcf ->
                snv: meta.variant_type == "snv"
                sv:  meta.variant_type == "sv"
            }
            .set { ch_mito_vcf_split }

        BCFTOOLS_VIEW_MITO.out.tbi
            .branch { meta, _tbi ->
                snv: meta.variant_type == "snv"
                sv:  meta.variant_type == "sv"
            }
            .set { ch_mito_tbi_split }

        ch_snv_vcf = ch_mito_vcf_split.snv.map { meta, vcf -> [meta - meta.subMap('variant_type'), vcf] }
        ch_snv_tbi = ch_mito_tbi_split.snv.map { meta, tbi -> [meta - meta.subMap('variant_type'), tbi] }
        ch_sv_vcf  = ch_mito_vcf_split.sv.map  { meta, vcf -> [meta - meta.subMap('variant_type'), vcf] }
        ch_sv_tbi  = ch_mito_tbi_split.sv.map  { meta, tbi -> [meta - meta.subMap('variant_type'), tbi] }

    } else {
        ch_snv_vcf = ch_vcf
        ch_snv_tbi = ch_tbi
        ch_sv_vcf  = channel.empty()
        ch_sv_tbi  = channel.empty()
    }

    emit:
    mitochondrial_snv_vcf = ch_snv_vcf  // channel: [val(meta), path(vcf)]
    mitochondrial_snv_tbi = ch_snv_tbi  // channel: [val(meta), path(tbi)]
    mitochondrial_sv_vcf  = ch_sv_vcf   // channel: [val(meta), path(vcf)]
    mitochondrial_sv_tbi  = ch_sv_tbi   // channel: [val(meta), path(tbi)]
}
