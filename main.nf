#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    genomic-medicine-sweden/nallo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/genomic-medicine-sweden/nallo
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NALLO                                 } from './workflows/nallo'
include { PIPELINE_INITIALISATION               } from './subworkflows/local/utils_nfcore_nallo_pipeline'
include { PIPELINE_COMPLETION                   } from './subworkflows/local/utils_nfcore_nallo_pipeline'
include {
    createReferenceChannelFromPath ;
    createReferenceChannelFromSamplesheet
} from './subworkflows/local/utils_nfcore_nallo_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow GENOMICMEDICINESWEDEN_NALLO {
    take:
    ch_cadd_header
    ch_cadd_prescored_indels
    ch_cadd_resources
    ch_echtvar_databases
    ch_exclude_bed
    ch_expected_xx_bed
    ch_expected_xy_bed
    ch_fasta
    ch_fai
    ch_genmod_reduced_penetrance
    ch_genmod_score_config_snvs
    ch_genmod_score_config_svs
    ch_gens_baf_positions
    ch_gens_coverage_bins
    ch_gens_panel_of_normals_female
    ch_gens_panel_of_normals_male
    ch_hgnc_ids
    ch_samplesheet
    ch_methbat_regions
    ch_modkit_call_regions
    ch_mosdepth_regions
    ch_paraphrase_rules
    ch_par
    ch_peddy_sites
    ch_sambamba_regions
    ch_sentieon_female_diploid_bed
    ch_sentieon_male_diploid_bed
    ch_sentieon_male_haploid_bed
    ch_sentieon_model_bundle
    ch_snv_call_regions
    ch_somalier_sites
    ch_stranger_repeat_catalog
    ch_str_bed
    ch_strdrop_training_set_json
    ch_sv_call_regions
    ch_svdb_sv_databases
    ch_tandem_repeats
    ch_variant_consequences_snvs
    ch_variant_consequences_svs
    ch_vcfexpress_prelude
    ch_vep_cache_unprocessed
    ch_vep_plugin_files
    cram_output
    val_alignment_processes
    val_bigwig_modcodes
    val_create_hificnv_maf_track
    val_create_sawfish_maf_track
    val_echtvar_snv_databases
    val_fasta
    val_filter_snvs_expression
    val_filter_svs_expression
    val_filter_variants_hgnc_ids
    val_force_sawfish_joint_call_single_samples
    val_hifiasm_mode
    val_mitochondrial_caller
    val_multiqc_config
    val_multiqc_logo
    val_multiqc_methods_description
    val_outdir
    val_paraphrase_output_format
    val_phaser
    val_plot_chromograph_autozygosity
    val_plot_chromograph_coverage
    val_pre_vep_snv_filter_expression
    val_run_methbat
    val_run_modkit
    val_sentieon_tech
    val_skip_alignment
    val_skip_annotate_paralogs
    val_skip_call_paralogs
    val_skip_chromograph
    val_skip_genome_assembly
    val_skip_methylation_calling
    val_skip_methylation_annotation
    val_skip_peddy
    val_skip_phasing
    val_skip_prepare_gens_input
    val_skip_qc
    val_skip_rank_variants
    val_skip_repeat_annotation
    val_skip_repeat_calling
    val_skip_sambamba_depth
    val_skip_sex_check
    val_skip_snv_annotation
    val_skip_snv_calling
    val_skip_sv_annotation
    val_skip_sv_calling
    val_snv_caller
    val_snv_calling_processes
    val_snv_call_regions
    val_str_caller
    val_strdrop_training_set_json
    val_sv_callers_merge_priority
    val_sv_callers_to_merge
    val_sv_callers_to_run
    val_sv_call_regions
    val_vep_cache
    val_vep_cache_version

    main:

    //
    // WORKFLOW: Run pipeline
    //
    NALLO(
        ch_cadd_header,
        ch_cadd_prescored_indels,
        ch_cadd_resources,
        ch_echtvar_databases,
        ch_exclude_bed,
        ch_expected_xx_bed,
        ch_expected_xy_bed,
        ch_fasta,
        ch_fai,
        ch_genmod_reduced_penetrance,
        ch_genmod_score_config_snvs,
        ch_genmod_score_config_svs,
        ch_gens_baf_positions,
        ch_gens_coverage_bins,
        ch_gens_panel_of_normals_female,
        ch_gens_panel_of_normals_male,
        ch_hgnc_ids,
        ch_samplesheet,
        ch_methbat_regions,
        ch_modkit_call_regions,
        ch_mosdepth_regions,
        ch_paraphrase_rules,
        ch_par,
        ch_peddy_sites,
        ch_sambamba_regions,
        ch_sentieon_female_diploid_bed,
        ch_sentieon_male_diploid_bed,
        ch_sentieon_male_haploid_bed,
        ch_sentieon_model_bundle,
        ch_snv_call_regions,
        ch_somalier_sites,
        ch_stranger_repeat_catalog,
        ch_str_bed,
        ch_strdrop_training_set_json,
        ch_sv_call_regions,
        ch_svdb_sv_databases,
        ch_tandem_repeats,
        ch_variant_consequences_snvs,
        ch_variant_consequences_svs,
        ch_vcfexpress_prelude,
        ch_vep_cache_unprocessed,
        ch_vep_plugin_files,
        cram_output,
        val_alignment_processes,
        val_bigwig_modcodes,
        val_create_hificnv_maf_track,
        val_create_sawfish_maf_track,
        val_echtvar_snv_databases,
        val_fasta,
        val_filter_snvs_expression,
        val_filter_svs_expression,
        val_filter_variants_hgnc_ids,
        val_force_sawfish_joint_call_single_samples,
        val_hifiasm_mode,
        val_mitochondrial_caller,
        val_multiqc_config,
        val_multiqc_logo,
        val_multiqc_methods_description,
        val_outdir,
        val_paraphrase_output_format,
        val_phaser,
        val_plot_chromograph_autozygosity,
        val_plot_chromograph_coverage,
        val_pre_vep_snv_filter_expression,
        val_run_methbat,
        val_run_modkit,
        val_sentieon_tech,
        val_skip_alignment,
        val_skip_annotate_paralogs,
        val_skip_call_paralogs,
        val_skip_chromograph,
        val_skip_genome_assembly,
        val_skip_methylation_calling,
        val_skip_methylation_annotation,
        val_skip_peddy,
        val_skip_phasing,
        val_skip_prepare_gens_input,
        val_skip_qc,
        val_skip_rank_variants,
        val_skip_repeat_annotation,
        val_skip_repeat_calling,
        val_skip_sambamba_depth,
        val_skip_sex_check,
        val_skip_snv_annotation,
        val_skip_snv_calling,
        val_skip_sv_annotation,
        val_skip_sv_calling,
        val_snv_caller,
        val_snv_calling_processes,
        val_snv_call_regions,
        val_str_caller,
        val_strdrop_training_set_json,
        val_sv_callers_merge_priority,
        val_sv_callers_to_merge,
        val_sv_callers_to_run,
        val_sv_call_regions,
        val_vep_cache,
        val_vep_cache_version,
    )

    emit:
    multiqc_report = NALLO.out.multiqc_report // channel: /path/to/multiqc_report.html
    aligned_assemblies = NALLO.out.aligned_assemblies // channel: [ val(meta), path(bam/cram), path(bai/crai) ]
    annotated_paralogs = NALLO.out.annotated_paralogs // channel: [ val(meta), path(tsv/json) ]
    annotated_repeats = NALLO.out.annotated_repeats // channel: [ val(meta), path(vcf), path(tbi) ]
    assembly_summary = NALLO.out.assembly_summary // channel: [ val(meta), path(assembly_summary) ]
    chromograph_plots = NALLO.out.chromograph_plots // channel: [ val(meta), path(png) ]
    gens_baf = NALLO.out.gens_baf // channel: [ val(meta), path(baf.bed.gz), path(baf.bed.gz.tbi) ]
    gens_cov = NALLO.out.gens_cov // channel: [ val(meta), path(cov.bed.gz), path(cov.bed.gz.tbi) ]
    haplotagged_reads = NALLO.out.haplotagged_reads // channel: [ val(meta), path(bam), path(bai) ]
    methylation_annotation = NALLO.out.methylation_annotation // channel: [ val(meta), path(methylated_regions_by_family) ]
    methylation_modkit_bed = NALLO.out.methylation_modkit_bed // channel: [ val(meta), path(bed.gz) ]
    methylation_modkit_tbi = NALLO.out.methylation_modkit_tbi // channel: [ val(meta), path(bed.gz.tbi) ]
    methylation_modkit_bigwig = NALLO.out.methylation_modkit_bigwig // channel: [ val(meta), path(bw) ]
    paralogs_family_vcf = NALLO.out.paralogs_family_vcf // channel: [ val(meta), path(vcf) ]
    paralogs_family_tbi = NALLO.out.paralogs_family_tbi // channel: [ val(meta), path(tbi) ]
    paralogs_sample_bam = NALLO.out.paralogs_sample_bam // channel: [ val(meta), path(bam) ]
    paralogs_sample_bai = NALLO.out.paralogs_sample_bai // channel: [ val(meta), path(bai) ]
    paralogs_sample_cram = NALLO.out.paralogs_sample_cram // channel: [ val(meta), path(cram) ]
    paralogs_sample_crai = NALLO.out.paralogs_sample_crai // channel: [ val(meta), path(crai) ]
    paralogs_sample_json = NALLO.out.paralogs_sample_json // channel: [ val(meta), path(json) ]
    paralogs_sample_vcf = NALLO.out.paralogs_sample_vcf // channel: [ val(meta), path(vcf) ]
    paralogs_sample_tbi = NALLO.out.paralogs_sample_tbi // channel: [ val(meta), path(tbi) ]
    repeat_trgt_sample_vcf = NALLO.out.repeat_trgt_sample_vcf // channel: [ val(meta), path(vcf) ]
    repeat_trgt_sample_tbi = NALLO.out.repeat_trgt_sample_tbi // channel: [ val(meta), path(tbi) ]
    repeat_trgt_sample_bam = NALLO.out.repeat_trgt_sample_bam // channel: [ val(meta), path(bam) ]
    repeat_trgt_sample_bai = NALLO.out.repeat_trgt_sample_bai // channel: [ val(meta), path(bai) ]
    repeat_trgt_sample_cram = NALLO.out.repeat_trgt_sample_cram // channel: [ val(meta), path(cram) ]
    repeat_trgt_sample_crai = NALLO.out.repeat_trgt_sample_crai // channel: [ val(meta), path(crai) ]
    repeat_trgt_family_vcf = NALLO.out.repeat_trgt_family_vcf // channel: [ val(meta), path(vcf) ]
    repeat_trgt_family_tbi = NALLO.out.repeat_trgt_family_tbi // channel: [ val(meta), path(tbi) ]
    cramino_unphased_stats = NALLO.out.cramino_unphased_stats // channel: [ val(meta), path(txt) ]
    cramino_unphased_arrow = NALLO.out.cramino_unphased_arrow // channel: [ val(meta), path(arrow) ]
    fastqc_html = NALLO.out.fastqc_html // channel: [ val(meta), path(html) ]
    fastqc_zip = NALLO.out.fastqc_zip // channel: [ val(meta), path(zip) ]
    mosdepth_summary = NALLO.out.mosdepth_summary // channel: [ val(meta), path(txt) ]
    mosdepth_global_dist = NALLO.out.mosdepth_global_dist // channel: [ val(meta), path(txt) ]
    mosdepth_regions_dist = NALLO.out.mosdepth_regions_dist // channel: [ val(meta), path(txt) ]
    mosdepth_per_base_d4 = NALLO.out.mosdepth_per_base_d4 // channel: [ val(meta), path(d4) ]
    mosdepth_regions = NALLO.out.mosdepth_regions // channel: [ val(meta), path(bed.gz), path(bed.gz.csi) ]
    sambamba_depth_bed = NALLO.out.sambamba_depth_bed // channel: [ val(meta), path(bed) ]
    qc_bcftools_stats = NALLO.out.qc_bcftools_stats // channel: [ val(meta), path(txt) ]
    qc_deepvariant_vcfstatsreport = NALLO.out.qc_deepvariant_vcfstatsreport // channel: [ val(meta), path(html) ]
    sample_snvs = NALLO.out.sample_snvs // channel: [ val(meta), path(vcf), path(tbi) ]
    somalier_relate_html    = NALLO.out.somalier_relate_html    // channel: [ val(meta), path(html) ]
    somalier_relate_pairs   = NALLO.out.somalier_relate_pairs   // channel: [ val(meta), path(pairs.tsv) ]
    somalier_relate_samples = NALLO.out.somalier_relate_samples // channel: [ val(meta), path(samples.tsv) ]
    phasing_stats = NALLO.out.phasing_stats // channel: [ val(meta), path("*.stats.tsv") ]
    phasing_blocks = NALLO.out.phasing_blocks // channel: [ val(meta), path("*.blocks.gtf.gz"), path("*.blocks.gtf.gz.tbi") ]
    haplotagging_stats = NALLO.out.haplotagging_stats // channel: [ val(meta), path("*.txt") ]
    haplotagging_arrow = NALLO.out.haplotagging_arrow // channel: [ val(meta), path("*.arrow") ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.help,
        params.help_full,
        args,
        params.outdir,
        params.show_hidden,
        params.validate_params,
        params.cnv_excluded_regions,
        params.cnv_expected_xx_cn,
        params.cnv_expected_xy_cn,
        params.echtvar_snv_databases,
        params.fasta,
        params.genmod_reduced_penetrance,
        params.genmod_score_config_snvs,
        params.genmod_score_config_svs,
        params.gens_baf_positions,
        params.gens_coverage_bins,
        params.gens_panel_of_normals_female,
        params.gens_panel_of_normals_male,
        params.input,
        params.methbat_regions,
        params.par_regions,
        params.phaser,
        params.run_methbat,
        params.sambamba_regions,
        params.skip_alignment,
        params.skip_annotate_paralogs,
        params.skip_call_paralogs,
        params.skip_chromograph,
        params.skip_genome_assembly,
        params.skip_methylation_calling,
        params.skip_methylation_annotation,
        params.skip_peddy,
        params.skip_phasing,
        params.skip_prepare_gens_input,
        params.skip_qc,
        params.skip_rank_variants,
        params.skip_repeat_annotation,
        params.skip_repeat_calling,
        params.skip_sambamba_depth,
        params.skip_sex_check,
        params.skip_snv_annotation,
        params.skip_snv_calling,
        params.skip_sv_annotation,
        params.skip_sv_calling,
        params.snv_caller,
        params.snv_calling_processes,
        params.somalier_sites,
        params.stranger_repeat_catalog,
        params.str_bed,
        params.str_caller,
        params.svdb_sv_databases,
        params.sv_callers,
        params.sv_callers_merge_priority,
        params.sv_callers_to_merge,
        params.sv_callers_to_run,
        params.variant_consequences_snvs,
        params.variant_consequences_svs,
        params.vep_cache,
        params.vep_plugin_files,
        params.version,
    )

    //
    // WORKFLOW: Run main workflow
    //
    GENOMICMEDICINESWEDEN_NALLO(
        createReferenceChannelFromPath("${projectDir}/assets/cadd_to_vcf_header_-1.0-.txt"),
        createReferenceChannelFromPath(params.cadd_prescored_indels),
        createReferenceChannelFromPath(params.cadd_resources),
        createReferenceChannelFromSamplesheet(params.echtvar_snv_databases, 'assets/schema_snp_db.json', channel.value([[], []])),
        createReferenceChannelFromPath(params.cnv_excluded_regions),
        createReferenceChannelFromPath(params.cnv_expected_xx_cn),
        createReferenceChannelFromPath(params.cnv_expected_xy_cn),
        createReferenceChannelFromPath(params.fasta),
        createReferenceChannelFromPath(params.fai),
        createReferenceChannelFromPath(params.genmod_reduced_penetrance),
        createReferenceChannelFromPath(params.genmod_score_config_snvs),
        createReferenceChannelFromPath(params.genmod_score_config_svs),
        createReferenceChannelFromPath(params.gens_baf_positions),
        createReferenceChannelFromPath(params.gens_coverage_bins),
        createReferenceChannelFromPath(params.gens_panel_of_normals_female, '', 'female_pon'),
        createReferenceChannelFromPath(params.gens_panel_of_normals_male, '', 'male_pon'),
        createReferenceChannelFromSamplesheet(params.filter_variants_hgnc_ids, 'assets/schema_hgnc_ids.json', channel.value([])).map { hgnc_id_list -> hgnc_id_list[0].toString() }.collectFile(name: 'hgnc_ids.txt', newLine: true, sort: true).map { file -> [[id: 'hgnc_ids'], file] }.collect(),
        PIPELINE_INITIALISATION.out.samplesheet,
        createReferenceChannelFromPath(params.methbat_regions),
        createReferenceChannelFromPath(params.modkit_call_regions, channel.value([[], []])),
        createReferenceChannelFromPath(params.mosdepth_regions, channel.value([[], []])),
        createReferenceChannelFromPath(params.paraphrase_rules, channel.value([[], []])),
        createReferenceChannelFromPath(params.par_regions),
        createReferenceChannelFromPath(params.peddy_sites, channel.value([[], []])),
        createReferenceChannelFromPath(params.sambamba_regions, channel.value([[], []])),
        createReferenceChannelFromPath(params.sentieon_female_diploid_bed, channel.value([[], []])),
        createReferenceChannelFromPath(params.sentieon_male_diploid_bed, channel.value([[], []])),
        createReferenceChannelFromPath(params.sentieon_male_haploid_bed, channel.value([[], []])),
        createReferenceChannelFromPath(params.sentieon_model_bundle, channel.value([[], []])),
        createReferenceChannelFromPath(params.snv_call_regions, channel.value([[], []])),
        createReferenceChannelFromPath(params.somalier_sites),
        createReferenceChannelFromPath(params.stranger_repeat_catalog),
        createReferenceChannelFromPath(params.str_bed),
        createReferenceChannelFromPath(params.strdrop_training_set_json),
        createReferenceChannelFromPath(params.sv_call_regions),
        createReferenceChannelFromSamplesheet(params.svdb_sv_databases, 'assets/svdb_query_vcf_schema.json', channel.value([])),
        createReferenceChannelFromPath(params.tandem_repeats, channel.value([[], []])),
        createReferenceChannelFromPath(params.variant_consequences_snvs),
        createReferenceChannelFromPath(params.variant_consequences_svs),
        file("${projectDir}/assets/vcf_express_found_in_prelude.lua"),
        createReferenceChannelFromPath(params.vep_cache, channel.value([[], []])),
        createReferenceChannelFromSamplesheet(params.vep_plugin_files, 'assets/schema_vep_plugin_files.json', channel.value([])),
        params.alignment_output_format == 'cram',
        params.alignment_processes,
        params.bigwig_modcodes,
        params.create_hificnv_maf_track,
        params.create_sawfish_maf_track,
        params.echtvar_snv_databases,
        params.fasta,
        params.filter_snvs_expression,
        params.filter_svs_expression,
        params.filter_variants_hgnc_ids,
        params.force_sawfish_joint_call_single_samples,
        params.hifiasm_mode,
        params.mitochondrial_caller,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
        params.paraphrase_output_format,
        params.phaser,
        params.plot_chromograph_autozygosity,
        params.plot_chromograph_coverage,
        params.pre_vep_snv_filter_expression,
        params.run_methbat,
        params.run_modkit,
        params.sentieon_tech,
        params.skip_alignment,
        params.skip_annotate_paralogs,
        params.skip_call_paralogs,
        params.skip_chromograph,
        params.skip_genome_assembly,
        params.skip_methylation_calling,
        params.skip_methylation_annotation,
        params.skip_peddy,
        params.skip_phasing,
        params.skip_prepare_gens_input,
        params.skip_qc,
        params.skip_rank_variants,
        params.skip_repeat_annotation,
        params.skip_repeat_calling,
        params.skip_sambamba_depth,
        params.skip_sex_check,
        params.skip_snv_annotation,
        params.skip_snv_calling,
        params.skip_sv_annotation,
        params.skip_sv_calling,
        params.snv_caller,
        params.snv_calling_processes,
        params.snv_call_regions,
        params.str_caller,
        params.strdrop_training_set_json,
        params.sv_callers_merge_priority,
        params.sv_callers_to_merge,
        params.sv_callers_to_run,
        params.sv_call_regions,
        params.vep_cache,
        params.vep_cache_version,
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        GENOMICMEDICINESWEDEN_NALLO.out.multiqc_report,
    )

    //
    // WORKFLOW OUTPUTS: Group files by publish directory
    //
    ch_gens = GENOMICMEDICINESWEDEN_NALLO.out.gens_baf
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.gens_cov)

    ch_methylation_modkit_pileup = GENOMICMEDICINESWEDEN_NALLO.out.methylation_modkit_bed
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.methylation_modkit_tbi)

    ch_qc_cramino_unphased = GENOMICMEDICINESWEDEN_NALLO.out.cramino_unphased_stats
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.cramino_unphased_arrow)

    ch_qc_cramino_phased = GENOMICMEDICINESWEDEN_NALLO.out.haplotagging_stats
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.haplotagging_arrow)

    ch_qc_phasing_stats = GENOMICMEDICINESWEDEN_NALLO.out.phasing_stats
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.phasing_blocks.map { meta, gtf, _tbi -> [meta, gtf] })
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.phasing_blocks.map { meta, _gtf, tbi -> [meta, tbi] })

    ch_qc_fastqc = GENOMICMEDICINESWEDEN_NALLO.out.fastqc_html
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.fastqc_zip)

    ch_qc_mosdepth = GENOMICMEDICINESWEDEN_NALLO.out.mosdepth_summary
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.mosdepth_global_dist)
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.mosdepth_regions_dist)
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.mosdepth_per_base_d4)
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.mosdepth_regions.map { meta, bed, _csi -> [meta, bed] })
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.mosdepth_regions.map { meta, _bed, csi -> [meta, csi] })

    ch_somalier_relate = GENOMICMEDICINESWEDEN_NALLO.out.somalier_relate_html
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.somalier_relate_samples)
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.somalier_relate_pairs)


    ch_paraphase_sample = GENOMICMEDICINESWEDEN_NALLO.out.paralogs_sample_json

    ch_paraphase_sample_bam = GENOMICMEDICINESWEDEN_NALLO.out.paralogs_sample_bam
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.paralogs_sample_bai)

    ch_paraphase_sample_cram = GENOMICMEDICINESWEDEN_NALLO.out.paralogs_sample_cram
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.paralogs_sample_crai)

    ch_paraphase_sample_vcfs = GENOMICMEDICINESWEDEN_NALLO.out.paralogs_sample_vcf
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.paralogs_sample_tbi)

    ch_paraphase_family = GENOMICMEDICINESWEDEN_NALLO.out.annotated_paralogs
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.paralogs_family_vcf)
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.paralogs_family_tbi)

    ch_repeats_sample_trgt = GENOMICMEDICINESWEDEN_NALLO.out.repeat_trgt_sample_vcf
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.repeat_trgt_sample_tbi)

    ch_repeats_sample_trgt_bam = GENOMICMEDICINESWEDEN_NALLO.out.repeat_trgt_sample_bam
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.repeat_trgt_sample_bai)

    ch_repeats_sample_trgt_cram = GENOMICMEDICINESWEDEN_NALLO.out.repeat_trgt_sample_cram
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.repeat_trgt_sample_crai)

    ch_repeats_family_trgt = GENOMICMEDICINESWEDEN_NALLO.out.repeat_trgt_family_vcf
        .mix(GENOMICMEDICINESWEDEN_NALLO.out.repeat_trgt_family_tbi)

    publish:
    aligned_assemblies = GENOMICMEDICINESWEDEN_NALLO.out.aligned_assemblies // channel: [ val(meta), path(bam/cram), path(bai/crai) ]
    paraphase_sample = ch_paraphase_sample // channel: [ val(meta), path(json) ]
    paraphase_sample_bam = ch_paraphase_sample_bam // channel: [ val(meta), path(bam/bai) ]
    paraphase_sample_cram = ch_paraphase_sample_cram // channel: [ val(meta), path(cram/crai) ]
    paraphase_sample_vcfs = ch_paraphase_sample_vcfs // channel: [ val(meta), path(vcf/tbi) ]
    paraphase_family = ch_paraphase_family // channel: [ val(meta), path(vcf/tbi/tsv/json) ]
    annotated_repeats = GENOMICMEDICINESWEDEN_NALLO.out.annotated_repeats // channel: [ val(meta), path(vcf), path(tbi) ]
    assembly_summary = GENOMICMEDICINESWEDEN_NALLO.out.assembly_summary // channel: [ val(meta), path(assembly_summary) ]
    chromograph_plots = GENOMICMEDICINESWEDEN_NALLO.out.chromograph_plots // channel: [ val(meta), path(png) ]
    gens = ch_gens // channel: [ val(meta), path(baf/cov.bed.gz), path(baf/cov.bed.gz.tbi) ]
    haplotagged_reads = GENOMICMEDICINESWEDEN_NALLO.out.haplotagged_reads // channel: [ val(meta), path(bam/cram), path(bai/crai) ]
    methylation_annotation = GENOMICMEDICINESWEDEN_NALLO.out.methylation_annotation // channel: [ val(meta), path(methylated_regions_by_family) ]
    methylation_modkit_pileup = ch_methylation_modkit_pileup // channel: [ val(meta), path(bed.gz/bed.gz.tbi) ]
    visualization_tracks_modkit = GENOMICMEDICINESWEDEN_NALLO.out.methylation_modkit_bigwig // channel: [ val(meta), path(bw) ]
    repeats_sample_trgt = ch_repeats_sample_trgt // channel: [ val(meta), path(vcf/tbi/bam/bai) ]
    repeats_sample_trgt_bam = ch_repeats_sample_trgt_bam  // channel: [ val(meta), path(bam/bai) ]
    repeats_sample_trgt_cram = ch_repeats_sample_trgt_cram // channel: [ val(meta), path(cram/crai) ]
    repeats_family_trgt = ch_repeats_family_trgt // channel: [ val(meta), path(vcf/tbi) ]
    qc_cramino_unphased = ch_qc_cramino_unphased // channel: [ val(meta), path(txt/arrow) ]
    qc_fastqc = ch_qc_fastqc // channel: [ val(meta), path(html/zip) ]
    qc_mosdepth = ch_qc_mosdepth // channel: [ val(meta), path(txt/d4/bed.gz/bed.gz.csi) ]
    qc_sambamba_depth = GENOMICMEDICINESWEDEN_NALLO.out.sambamba_depth_bed // channel: [ val(meta), path(bed) ]
    qc_bcftools_stats = GENOMICMEDICINESWEDEN_NALLO.out.qc_bcftools_stats // channel: [ val(meta), path(txt) ]
    qc_deepvariant_vcfstatsreport = GENOMICMEDICINESWEDEN_NALLO.out.qc_deepvariant_vcfstatsreport // channel: [ val(meta), path(html) ]
    sample_snvs = GENOMICMEDICINESWEDEN_NALLO.out.sample_snvs // channel: [ val(meta), path(vcf), path(tbi) ]
    somalier_relate = ch_somalier_relate // channel: [ val(meta), path(html/pairs/samples) ]
    qc_cramino_phased = ch_qc_cramino_phased // channel: [ val(meta), path(txt/arrow) ]
    qc_phasing_stats = ch_qc_phasing_stats // channel: [ val(meta), path(tsv/gtf.gz/gtf.gz.tbi) ]
}

output {
    aligned_assemblies {
        path { meta, _bam, _bai -> "assembly/sample/${meta.id}/" }
    }
    paraphase_sample {
        path { meta, _file -> "paraphase/sample/${meta.id}/" }
    }
    paraphase_sample_bam {
        path { meta, _file -> "paraphase/sample/${meta.id}/" }
        enabled params.alignment_output_format == 'bam'
    }
    paraphase_sample_cram {
        path { meta, _file -> "paraphase/sample/${meta.id}/" }
        enabled params.alignment_output_format == 'cram'
    }
    paraphase_sample_vcfs {
        path { meta, _file -> "paraphase/sample/${meta.id}/" }
    }
    paraphase_family {
        path { meta, _file -> "paraphase/family/${meta.id}/" }
    }
    annotated_repeats {
        path { meta, _vcf, _tbi -> "repeats/family/${meta.id}/" }
    }
    assembly_summary {
        path { meta, _assembly_summary -> "assembly/stats/${meta.id}/" }
    }
    chromograph_plots {
        path { meta, _plot -> "images/chromograph/sample/${meta.id}/" }
    }
    gens {
        path { meta, _bed, _tbi -> "gens/${meta.id}/" }
    }
    haplotagged_reads { // HiPhase uses the input file (aligned reads) as template for naming output, so we need to remove the "_aligned" suffix here
        path { meta, bam, bai ->
            bam >> "aligned_reads/${meta.id}/${bam.name.replaceFirst("_aligned", "")}"
            bai >> "aligned_reads/${meta.id}/${bai.name.replaceFirst("_aligned", "")}"
            }
    }
    methylation_annotation {
        path { meta, _methylated_regions -> "methylation/profile/family/${meta.id}/" }
    }
    methylation_modkit_pileup {
        path { meta, _file -> "methylation/pileup/${meta.id}/" }
    }
    visualization_tracks_modkit {
        path { meta, _bw -> "visualization_tracks/${meta.id}/" }
    }
    repeats_sample_trgt {
        path { meta, _file -> "repeats/sample/${meta.id}/" }
    }
    repeats_sample_trgt_cram {
        path { meta, _file -> "repeats/sample/${meta.id}/" }
        enabled params.alignment_output_format == 'cram'
    }
    repeats_sample_trgt_bam {
        path { meta, _file -> "repeats/sample/${meta.id}/" }
        enabled params.alignment_output_format == 'bam'
    }
    repeats_family_trgt {
        path { meta, _file -> "repeats/family/${meta.id}/" }
    }
    qc_cramino_unphased {
        path { meta, _file -> "qc/cramino/unphased/${meta.id}/" }
    }
    qc_fastqc {
        path { meta, _file -> "qc/fastqc/${meta.id}/" }
    }
    qc_mosdepth {
        path { meta, _file -> "qc/mosdepth/${meta.id}/" }
    }
    qc_sambamba_depth {
        path { meta, _file -> "qc/sambamba_depth/${meta.id}/" }
    }
    qc_bcftools_stats {
        path { meta, _stats -> "qc/bcftools_stats/${meta.id}/" }
    }
    qc_deepvariant_vcfstatsreport {
            path { meta, _report -> "qc/deepvariant_vcfstatsreport/${meta.id}/" }
    }
    qc_cramino_phased {
        path { meta, _file -> "qc/cramino/phased/${meta.id}/" }
    }
    qc_phasing_stats {
        path { meta, _file -> "qc/phasing_stats/${meta.id}/" }
    }
    sample_snvs {
        path { meta, _vcf, _tbi -> "snvs/sample/${meta.id}/" }
    }
    somalier_relate {
        path { meta, _file -> "qc/somalier/relate/${meta.id}/" }
    }
}
