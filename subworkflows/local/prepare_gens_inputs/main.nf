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
    ch_gvcfs            // channel: [mandatory] [ val(meta), path(gvcf) ]
    baf_positions       // channel: [mandatory] [ path(gz) ]
    gatk_header         // channel: [mandatory] [ path(txt) ]
    panel_of_normals    // channel: [mandatory] [ path(hd5) ]

    main:
    ch_versions = channel.empty()

    //ch_gvcfs.view()

    // ChatGPT in progress
    // 1. Flatten ?
    // 2. Group ?

    ch_remapped = ch_gvcfs.map { meta, pairs -> 
        def vcfs = pairs.collect { it -> it[0] }
        def tbis = pairs.collect { it -> it[1] }
        tuple(meta, vcfs, tbis)
    }


    // ch_remapped.view { it -> "ch_remapped ${it}" }

    // FIXME: Looks like the nesting is not correct yet [meta, (vcf, tbi), (vcf, tbi)] vs [meta, (vcfs), (tbi)]
    BCFTOOLS_CONCAT(ch_remapped)

    ch_vcf_tbi = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi)
    // ch_vcf_tbi.view()

    // BCFTOOLS_CONCAT.out.vcf.view()
    // BCFTOOLS_CONCAT.out.tbi.view()

    ch_mosdepth_in = ch_bam.map { meta, bam, bai -> tuple(meta, bam, bai, []) }
    ch_empty_fasta = ch_bam.map { meta, _bam, _bai -> tuple(meta, []) }

    MOSDEPTH(ch_mosdepth_in, ch_empty_fasta)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    MOSDEPTH_TO_GATK_FORMAT(
        MOSDEPTH.out.regions_bed,
        gatk_header,
    )
    ch_versions = ch_versions.mix(MOSDEPTH_TO_GATK_FORMAT.out.versions)

    MOSDEPTH_TO_GATK_FORMAT.out.output
        .map { meta, _tsv -> tuple(meta, panel_of_normals) }
        .set { ch_pon }

    //MOSDEPTH_TO_GATK_FORMAT.out.output
    //    .map { meta, _tsv -> meta }
    //    .cross(panel_of_normals)
    //    .map { meta, pon -> tuple(meta, pon) }
    //    .set { ch_pon }

    GATK4_DENOISEREADCOUNTS(
        MOSDEPTH_TO_GATK_FORMAT.out.output,
        ch_pon
        // meta, pon
    )
    ch_versions = ch_versions.mix(GATK4_DENOISEREADCOUNTS.out.versions)

    ch_gens_input = GATK4_DENOISEREADCOUNTS.out.standardized
        .join(ch_vcf_tbi)

    ch_gens_input.view()

    // Input: meta, coverage, gvcf
    GENERATE_GENS_DATA(
        ch_gens_input,
        baf_positions
    )
    ch_versions = ch_versions.mix(GENERATE_GENS_DATA.out.versions)

    // FIXME: Do we want a "without normal" cov output as well?
    emit:
    cov_bed_tbi = GENERATE_GENS_DATA.out.cov_bed_tbi
    baf_bed_tbi = GENERATE_GENS_DATA.out.baf_bed_tbi
    versions = ch_versions
}
