include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'
include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'
include { GENERATE_GENS_DATA } from '../../../modules/local/generate_gens_data/main'
include { MOSDEPTH_TO_GATK_FORMAT } from '../../../modules/local/mosdepth_to_gatk_format/main'
include { GATK4_DENOISEREADCOUNTS } from '../../../modules/nf-core/gatk4/denoisereadcounts/main'

// Calculate the GATK header from this instead
// samtools view -H sample.bam > header.sam


// FIXME: Dealing with multiple samples
workflow PREPARE_GENS_INPUTS {
    take:
    ch_bam              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_gvcfs            // channel: [mandatory] [ val(meta), tuple(path(gvcfs)), tuple(path(tbis))]
    ch_baf_positions       // channel: [mandatory] [ path(gz) ]
    // FIXME: Calculate from BAM instead
    gatk_header         // value: [mandatory] [ path(txt) ]
    ch_panel_of_normals    // channel: [mandatory] [ path(hd5) ]

    main:
    ch_versions = channel.empty()

    // These looks OK
    //ch_bam.view()
    //ch_gvcfs.view()
    //ch_baf_positions.view()
    //ch_gatk_header.view()
    //ch_panel_of_normals.view()

    BCFTOOLS_CONCAT(ch_gvcfs)
    ch_vcf_tbi = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi)

    ch_mosdepth_in = ch_bam.map { meta, bam, bai -> tuple(meta, bam, bai, []) }
    ch_empty_fasta = ch_bam.map { meta, _bam, _bai -> tuple(meta, []) }

    MOSDEPTH(ch_mosdepth_in, ch_empty_fasta)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    MOSDEPTH.out.regions_bed.view { it -> "Mosdepth: ${it}" }
    //ch_gatk_header.view { it -> "gatk header: ${it}" }

    MOSDEPTH.out.regions_bed
        .map { meta, regions_bed -> tuple(meta, regions_bed, gatk_header) }
        .set { ch_mosdepth_to_gatk_in }

    ch_mosdepth_to_gatk_in.view { it -> "ch_mosdepth_to_gatk_in ${it}" }

    MOSDEPTH_TO_GATK_FORMAT(ch_mosdepth_to_gatk_in)
    ch_versions = ch_versions.mix(MOSDEPTH_TO_GATK_FORMAT.out.versions)

    MOSDEPTH_TO_GATK_FORMAT.out.output.view()

    //MOSDEPTH_TO_GATK_FORMAT.out.output
    //    .combine(ch_panel_of_normals)
    //    .map { meta, _tsv, pon -> tuple(meta, pon) }
    //    .set { ch_pon }

    //MOSDEPTH_TO_GATK_FORMAT.out.output
    //    .map { meta, _tsv -> meta }
    //    .cross(panel_of_normals)
    //    .map { meta, pon -> tuple(meta, pon) }
    //    .set { ch_pon }

    //MOSDEPTH_TO_GATK_FORMAT.out.output.view()
    //ch_pon.view()
    //GATK4_DENOISEREADCOUNTS(
    //    MOSDEPTH_TO_GATK_FORMAT.out.output,
    //    ch_pon
    //    // meta, pon
    //)
    //ch_versions = ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions)

    //GATK4_DENOISEREADCOUNTS.out.standardized.view()

    //ch_gens_input = GATK4_DENOISEREADCOUNTS.out.standardized
    //    .join(ch_vcf_tbi)

    //ch_gens_input.view()

    // Input: meta, coverage, gvcf
    //GENERATE_GENS_DATA(
    //    ch_gens_input,
    //    ch_baf_positions
    //)
    //ch_versions = ch_versions.mix(GENERATE_GENS_DATA.out.versions)

    // FIXME: Do we want a "without normal" cov output as well?
    emit:
    //cov_bed_tbi = GENERATE_GENS_DATA.out.cov_bed_tbi
    //baf_bed_tbi = GENERATE_GENS_DATA.out.baf_bed_tbi
    versions = ch_versions
}
