include { samplesheetToList                                      } from 'plugin/nf-schema'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGN_ASSEMBLIES                                       } from '../subworkflows/local/align_assemblies'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV                    } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SVS                    } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_PARALOGS                                      } from '../subworkflows/local/annotate_paralogs'
include { ANNOTATE_REPEAT_EXPANSIONS                             } from '../subworkflows/local/annotate_repeat_expansions'
include { ANNOTATE_SNVS                                          } from '../subworkflows/local/annotate_snvs'
include { ANNOTATE_SVS                                           } from '../subworkflows/local/annotate_svs'
include { CONVERT_INPUT_FILES as CONVERT_INPUT_FASTQS            } from '../subworkflows/local/convert_input_files'
include { CONVERT_INPUT_FILES as CONVERT_INPUT_BAMS              } from '../subworkflows/local/convert_input_files'
include { CHROMOGRAPH                                            } from '../subworkflows/local/chromograph'
include { BAM_INFER_SEX                                          } from '../subworkflows/local/bam_infer_sex'
include { CALL_PARALOGS                                          } from '../subworkflows/local/call_paralogs'
include { CALL_REPEAT_EXPANSIONS_STRDUST                         } from '../subworkflows/local/call_repeat_expansions_strdust'
include { CALL_REPEAT_EXPANSIONS_TRGT                            } from '../subworkflows/local/call_repeat_expansions_trgt'
include { CALL_SNVS                                              } from '../subworkflows/local/call_snvs'
include { CALL_SVS                                               } from '../subworkflows/local/call_svs'
include { GENOME_ASSEMBLY                                        } from '../subworkflows/local/genome_assembly'
include { GVCF_GLNEXUS_NORM_VARIANTS                             } from '../subworkflows/local/gvcf_glnexus_norm_variants'
include { CALL_METHYLATION_MODKIT                                } from '../subworkflows/local/call_methylation_modkit'
include { CALL_METHYLATION_METHBAT                               } from '../subworkflows/local/call_methylation_methbat'
include { PHASING                                                } from '../subworkflows/local/phasing'
include { PREPARE_GENS_INPUTS                                    } from '../subworkflows/local/prepare_gens_inputs'
include { PREPARE_REFERENCES                                     } from '../subworkflows/local/prepare_references'
include { QC_ALIGNED_READS                                       } from '../subworkflows/local/qc_aligned_reads'
include { QC_SNVS                                                } from '../subworkflows/local/qc_snvs'
include { RANK_VARIANTS as RANK_VARIANTS_SNV                     } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SVS                     } from '../subworkflows/local/rank_variants'
include { SCATTER_GENOME                                         } from '../subworkflows/local/scatter_genome'
include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP as FILTER_VARIANTS_SNVS } from '../subworkflows/nf-core/vcf_filter_bcftools_ensemblvep/main'
include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP as FILTER_VARIANTS_SVS  } from '../subworkflows/nf-core/vcf_filter_bcftools_ensemblvep/main'
include { VCF_CONCAT_NORM_VARIANTS                               } from '../subworkflows/local/vcf_concat_norm_variants'
include { VCF_CONCAT_SORT_VARIANTS as CONCAT_SORT_ANNOTATED_SNVS } from '../subworkflows/local/vcf_concat_sort_variants/main'
include { VCF_CONCAT_SORT_VARIANTS as CONCAT_SORT_RANKED_SNVS    } from '../subworkflows/local/vcf_concat_sort_variants/main'
include { VCF_CONCAT_SORT_VARIANTS as CONCAT_SORT_GENS           } from '../subworkflows/local/vcf_concat_sort_variants/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// local
include { CREATE_PEDIGREE_FILE as SAMPLESHEET_PED                } from '../modules/local/create_pedigree_file/main'
include { CREATE_PEDIGREE_FILE as SOMALIER_PED_FAMILY            } from '../modules/local/create_pedigree_file/main'

// nf-core
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_PHASING             } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CHROMOGRAPH             } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SV                      } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_PHASING                 } from '../modules/nf-core/bcftools/view/main'
include { MINIMAP2_ALIGN                                         } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                                         } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_CONVERT                                       } from '../modules/nf-core/samtools/convert/main'
include { MULTIQC                                                } from '../modules/nf-core/multiqc/main'
include { PEDDY                                                  } from '../modules/nf-core/peddy/main'
include { SPLITUBAM                                              } from '../modules/nf-core/splitubam/main'
include { SVDB_MERGE as SVDB_MERGE_SVS_CNVS                      } from '../modules/nf-core/svdb/merge/main'
include { paramsSummaryMap                                       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                 } from '../subworkflows/local/utils_nfcore_nallo_pipeline'
include { citationBibliographyText                               } from '../subworkflows/local/utils_nfcore_nallo_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NALLO {
    take:
    ch_cadd_header
    ch_cadd_prescored_indels
    ch_cadd_resources
    ch_echtvar_databases
    ch_exclude_bed
    ch_expected_xx_bed
    ch_expected_xy_bed
    ch_fasta
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
    ch_multiqc_files = channel.empty()


    //
    // Prepare references
    //
    if (!val_skip_alignment || !val_skip_genome_assembly) {

        // The genome assembly alignment needs a fai for cram output,
        // but we shouldn't need to prepare the vep cache here.
        // Perhaps PREPARE_REFERENCES could be modified to handle this case?
        PREPARE_REFERENCES(
            ch_fasta,
            ch_vep_cache_unprocessed,
            val_fasta.endsWith('.gz'),
            val_vep_cache && val_vep_cache.endsWith("tar.gz"),
        )

        // Gather indices
        ch_fasta = PREPARE_REFERENCES.out.fasta
        ch_fai = PREPARE_REFERENCES.out.fai
    }

    // Convert FASTQ to BAM only if alignment or should be done.
    // Since we assume that the majority of pipeline runs will use BAM files as input,
    // we start all files as BAMs for simplicity except for the assembly, which requires FASTQs.
    if (!val_skip_alignment) {

        CONVERT_INPUT_FASTQS(
            ch_samplesheet,
            false,
            true,
        )
    }

    // To speed up the alignement, we can split the BAM files into smaller chunks.
    // We can also use the split BAM files for FASTQ conversion to the assembly workflow,
    // instead of the original BAM files which should allow the assembly to start sooner.
    //
    // We could change the name of alignment processes to something more generic, like `--split_input_files`?
    // If we start running more trios we also need to consider that the parents at the moment needs to be merged
    // before YAK. So, we could consider adding some logic to handle that case,
    // to avoid unneccessary splitting and merging just for a minor speedup in the conversion.
    if (!val_skip_alignment && val_alignment_processes > 1) {

        // contains all BAM files, including those not converted.
        SPLITUBAM(
            CONVERT_INPUT_FASTQS.out.bam
        )
    }

    //
    // Hifiasm assembly and alignment to reference genome
    //
    if (!val_skip_genome_assembly) {

        // Now, if we started with BAM files, we do alignment and split the BAM files (this is the main case),
        // then we converted any FASTQs to BAMs, the original BAMs will have been split, and we can convert those
        // to FASTQs for the assembly.
        //
        // We could possibly implement something where we would check if the converted BAM had an original FASTQ,
        // then we could use that FASTQ directly for the assembly. But since this is a rare case, it's not implemented.
        //
        // If we didn't split the files, there's currently no need to take the converted BAMs,
        // so we take the original FASTQs instead if there are any, and also convert the original BAMs to FASTQs,
        // so we can use those for the assembly.
        //
        // Since starting with FASTQs is a rare case, no splitting of FASTQs alone just for the assembly is implmenented

        CONVERT_INPUT_BAMS(
            !val_skip_alignment && val_alignment_processes > 1 ? SPLITUBAM.out.bam.transpose() : ch_samplesheet,
            true,
            false,
        )

        // contains all FASTQ files, including those not converted
        CONVERT_INPUT_BAMS.out.fastq
            .groupTuple()
            .set { ch_genome_assembly_input }

        // Hifiasm assembly
        GENOME_ASSEMBLY(
            ch_genome_assembly_input,
            val_hifiasm_mode == "trio-binning",
        )

        ALIGN_ASSEMBLIES(
            GENOME_ASSEMBLY.out.assembled_haplotypes,
            ch_fasta,
            ch_fai,
            cram_output,
        )
    }

    /*
     * Map reads to reference
     */
    if (!val_skip_alignment) {

        /*
         * Ensure each BAM has a unique identify,
         * enabling correct grouping and downstream merging.
         */
        (val_alignment_processes > 1
            ? SPLITUBAM.out.bam.transpose()
            : CONVERT_INPUT_FASTQS.out.bam).map { meta, bam -> [meta + [file: bam.name], bam] }.set { reads_for_alignment }

        /*
         * Create a grouping key per sample that records the number of split files,
         * allowing downstream merging to trigger as soon as all alignments of a sample are ready.
         */
        reads_for_alignment
            .map { meta, bam -> [meta - meta.subMap('file'), bam] }
            .groupTuple()
            .map { meta, files -> [meta + [n_files: files.size()]] }
            .set { reads_grouping_key }

        /*
         * Align reads independently per split (could be a split-align-merge subworkflow)
         */
        MINIMAP2_ALIGN(
            reads_for_alignment,
            PREPARE_REFERENCES.out.mmi,
            true,
            'bai',
            false,
            false,
        )

        /*
         * Re-attach grouping key so BAMs can be merged per group as soon as all alignments for one sample are ready
         */
        MINIMAP2_ALIGN.out.bam
            .join(MINIMAP2_ALIGN.out.index, failOnDuplicate: true, failOnMismatch: true)
            .combine(reads_grouping_key)
            .filter { bam_meta, _bam, _bai, grouping_key_meta ->
                bam_meta.id == grouping_key_meta.id
            }
            .map { bam_meta, bam, bai, grouping_key_meta ->
                [bam_meta - bam_meta.subMap('file') + [n_files: grouping_key_meta.n_files], bam, bai]
            }
            .map { meta, bam, bai ->
                [groupKey(meta, meta.n_files), bam, bai]
            }
            .groupTuple()
            .set { bam_to_merge }

        /*
         * Always merge here even if there's only one file,
         * because alignment runs without knowledge of group completeness (n_files),
         * and we can't therefore output from the alignment step with correct naming.
         */
        SAMTOOLS_MERGE(
            bam_to_merge,
            [[], [], [], []],
        )

        // Combine merged with unmerged bam files
        SAMTOOLS_MERGE.out.bam
            .join(SAMTOOLS_MERGE.out.index, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, bam, bai -> [meta - meta.subMap('n_files'), bam, bai] }
            .set { ch_aligned_bam }

        // Publish alignments as CRAM if requested
        if (cram_output && val_skip_phasing) {
            SAMTOOLS_CONVERT(
                ch_aligned_bam,
                ch_fasta.join(ch_fai).collect(),
            )
        }

        //
        // Create PED from samplesheet
        //
        ch_samplesheet
            .map { meta, _files -> [[id: meta.project], meta] }
            .groupTuple()
            .set { ch_samplesheet_ped_in }

        SAMPLESHEET_PED(ch_samplesheet_ped_in)

        SAMPLESHEET_PED.out.ped
            .collect()
            .set { ch_samplesheet_pedfile }

        if (!val_skip_sex_check) {
            //
            // Check sex and relatedness, and update with inferred sex if the sex for a sample is unknown
            //
            BAM_INFER_SEX(
                ch_aligned_bam,
                ch_fasta,
                ch_fai,
                ch_somalier_sites,
                ch_samplesheet_pedfile,
            )
            ch_multiqc_files = ch_multiqc_files.mix(BAM_INFER_SEX.out.somalier_samples.map { _meta, metrics -> metrics }.collect().ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(BAM_INFER_SEX.out.somalier_pairs.map { _meta, metrics -> metrics }.collect().ifEmpty([]))

            // Set files with updated meta for subsequent processes
            ch_bam = BAM_INFER_SEX.out.bam
            ch_bam_bai = BAM_INFER_SEX.out.bam_bai
        }
        else {
            ch_bam = ch_aligned_bam.map { meta, bam, _bai -> [meta, bam] }
            ch_bam_bai = ch_aligned_bam
        }
    }

    //
    // Run read QC with FastQC, mosdepth and cramino
    //
    if (!val_skip_qc) {

        QC_ALIGNED_READS(
            ch_bam_bai,
            ch_fasta,
            ch_mosdepth_regions,
            ch_sambamba_regions,
            !val_skip_sambamba_depth,
        )
        ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNED_READS.out.fastqc_zip.collect { _meta, metrics -> metrics }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNED_READS.out.mosdepth_summary.collect { _meta, metrics -> metrics })
        ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNED_READS.out.mosdepth_global_dist.collect { _meta, metrics -> metrics })
        ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNED_READS.out.mosdepth_region_dist.collect { _meta, metrics -> metrics }.ifEmpty([]))
    }

    /*
     * Call paralogous genes with paraphase
     */
    if (!val_skip_call_paralogs) {
        CALL_PARALOGS(
            ch_bam_bai,
            ch_fasta,
            ch_fai,
            cram_output,
        )
    }

    /*
     * Annotate paralogous genes with paraphrase
     */
    if (!val_skip_annotate_paralogs) {
        ANNOTATE_PARALOGS(
            CALL_PARALOGS.out.json,
            val_paraphrase_output_format,
            ch_paraphrase_rules,
        )
    }

    /*
     * Call SNVs
     */
    if (!val_skip_snv_calling) {

        // Make BED intervals, can be used for parallel SNV calling
        SCATTER_GENOME(
            ch_fai,
            ch_snv_call_regions,
            !val_snv_call_regions,
            val_snv_calling_processes,
        )

        if (val_mitochondrial_caller == "deepvariant") {

            SCATTER_GENOME.out.bed_nuclear_mitochondrial_intervals.set { ch_bed_intervals }
        }
        else {
            ch_bed_intervals = SCATTER_GENOME.out.bed_nuclear_intervals
        }

        // Combine the BED intervals with BAM/BAI files to create a region-bam-bai for each sample.
        // This uses the whole BAM files for each region instead of splitting them.
        ch_bam_bai
            .combine(ch_bed_intervals)
            .map { meta, bam, bai, bed_meta, bed, num_intervals ->
                [meta + [genome: bed_meta.genome, num_intervals: num_intervals, region: bed], bam, bai, bed]
            }
            .set { call_snvs_input }

        CALL_SNVS(
            call_snvs_input,
            ch_fasta,
            ch_fai,
            ch_par,
            ch_sentieon_model_bundle,
            ch_sentieon_female_diploid_bed,
            ch_sentieon_male_diploid_bed,
            ch_sentieon_male_haploid_bed,
            val_snv_caller,
            val_sentieon_tech,
        )

        // Group GVCFs per region and family (one region with all samples)
        CALL_SNVS.out.gvcf
            .map { meta, gvcf ->
                [[id: meta.region.name, family_id: meta.family_id, num_intervals: meta.num_intervals], gvcf]
            }
            .groupTuple()
            .set { variants_to_merge_per_family }

        CALL_SNVS.out.gvcf_index
            .map { meta, tbi ->
                [[id: meta.region.name, family_id: meta.family_id, num_intervals: meta.num_intervals], tbi]
            }
            .groupTuple()
            .set { gvcf_tbis_per_family }

        // Create a merged and normalized VCF, containing one region with all samples, to be used in annotation and ranking.
        // SCATTER_GENOME.out.bed contains all regions, but we could probably pass the region BED that actually matches the variants instead...
        GVCF_GLNEXUS_NORM_VARIANTS(
            variants_to_merge_per_family,
            gvcf_tbis_per_family,
            SCATTER_GENOME.out.bed,
            ch_fasta,
            ch_fai,
            val_snv_caller,
            ch_vcfexpress_prelude,
        )

        // Grouping VCF, containing one sample with all regions
        CALL_SNVS.out.vcf
            .map { meta, vcf ->
                def new_meta = meta - meta.subMap('region', 'genome')
                [groupKey(new_meta, new_meta.num_intervals), vcf]
            }
            .groupTuple()
            .map { meta, vcfs ->
                [meta - meta.subMap('num_intervals'), vcfs]
            }
            .set { variants_to_concat_per_sample }

        // Create a concatenated and normalized VCF, containing one sample with all regions.
        VCF_CONCAT_NORM_VARIANTS(
            variants_to_concat_per_sample,
            ch_fasta,
            val_snv_caller,
            ch_vcfexpress_prelude,
        )

        // These contains RefCalls
        sample_snv_vcf = VCF_CONCAT_NORM_VARIANTS.out.vcf
        sample_snv_index = VCF_CONCAT_NORM_VARIANTS.out.index

        family_snv_vcf = GVCF_GLNEXUS_NORM_VARIANTS.out.vcf
        family_snv_index = GVCF_GLNEXUS_NORM_VARIANTS.out.index

        // SNV QC
        // Can we use the normalized VCF here, for DV vcfstatsreport?
        QC_SNVS(
            VCF_CONCAT_NORM_VARIANTS.out.bcftools_concat_vcf,
            sample_snv_vcf,
            sample_snv_index,
            val_snv_caller.equals("deepvariant"),
        )
        ch_multiqc_files = ch_multiqc_files.mix(QC_SNVS.out.stats.collect { _meta, metrics -> metrics }.ifEmpty([]))

        family_snv_vcf
            .join(family_snv_index, failOnMismatch: true, failOnDuplicate: true)
            .set { ch_vcf_tbi_per_region }
    }
    if (!val_skip_prepare_gens_input) {
        CALL_SNVS.out.gvcf
            .join(CALL_SNVS.out.gvcf_index)
            .map { meta, gvcf, gvcf_index ->
                def sample_meta = meta - meta.subMap(['region', 'num_intervals', 'genome'])
                [sample_meta, gvcf, gvcf_index]
            }
            .groupTuple()
            .set { ch_gvcfs }

        CONCAT_SORT_GENS(
            ch_gvcfs
        )

        CONCAT_SORT_GENS.out.vcf
            .join(CONCAT_SORT_GENS.out.index)
            .set { ch_gvcf_tbi }

        PREPARE_GENS_INPUTS(
            ch_bam_bai,
            ch_gvcf_tbi,
            ch_gens_baf_positions,
            ch_gens_panel_of_normals_female,
            ch_gens_panel_of_normals_male,
            ch_gens_coverage_bins,
        )
    }

    //
    // Call SVs
    //
    if (!val_skip_sv_calling) {

        CALL_SVS(
            ch_bam_bai,
            ch_tandem_repeats,
            sample_snv_vcf,
            ch_fasta,
            ch_expected_xy_bed,
            ch_expected_xx_bed,
            ch_exclude_bed,
            val_sv_callers_to_run.split(',').collect { caller -> caller.toLowerCase().trim() },
            val_sv_callers_to_merge.split(',').collect { caller -> caller.toLowerCase().trim() },
            val_sv_callers_merge_priority.split(',').collect { caller -> caller.toLowerCase().trim() },
            ch_sv_call_regions,
            val_sv_call_regions,
            val_force_sawfish_joint_call_single_samples,
            val_create_hificnv_maf_track,
            val_create_sawfish_maf_track,
            ch_vcfexpress_prelude,
        )
    }

    //
    // Phase SNVs and INDELs
    //
    if (!val_skip_phasing) {

        ch_samplesheet
            .map { meta, _files -> [meta.family_id, meta.id] }
            .groupTuple()
            .map { family_id, sample_ids ->
                [[id: family_id], sample_ids.unique()]
            }
            .set { ch_family_to_samples }

        // Grouping SNV VCFs per family to concatenate before phasing.
        // Right now they are split by calling regions but we need whole-genome VCFs for phasing.
        family_snv_vcf
            .join(family_snv_index, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, vcf, tbi ->
                def new_meta = [id: meta.family_id, num_intervals: meta.num_intervals]
                [groupKey(new_meta, meta.num_intervals), vcf, tbi]
            }
            .groupTuple()
            .map { key, vcfs, tbis ->
                [key.getGroupTarget(), vcfs, tbis]
            }
            .map { meta, vcfs, tbis ->
                [meta - meta.subMap('num_intervals'), vcfs, tbis]
            }
            .set { ch_bcftools_concat_phasing_in }

        BCFTOOLS_CONCAT_PHASING(
            ch_bcftools_concat_phasing_in
        )

        PHASING(
            BCFTOOLS_CONCAT_PHASING.out.vcf,
            BCFTOOLS_CONCAT_PHASING.out.tbi,
            val_skip_sv_calling ? channel.empty() : CALL_SVS.out.family_vcf,
            val_skip_sv_calling ? channel.empty() : CALL_SVS.out.family_tbi,
            ch_bam_bai,
            ch_family_to_samples,
            ch_fasta,
            ch_fai,
            val_phaser,
            !val_skip_sv_calling,
            cram_output,
        )

        ch_multiqc_files = ch_multiqc_files.mix(PHASING.out.stats.collect { _meta, txt -> txt }.ifEmpty([]))

        // Scatter whole-genome phased SNV VCFs back into regions for annotation
        PHASING.out.phased_family_snvs
            .join(PHASING.out.phased_family_snvs_tbi, failOnMismatch: true, failOnDuplicate: true)
            .combine(ch_bed_intervals)
            .multiMap { vcf_meta, vcf, tbi, bed_meta, bed, num_intervals ->
                vcf: [vcf_meta + [id: bed.name, family_id: vcf_meta.id, genome: bed_meta.genome, num_intervals: num_intervals], vcf, tbi]
                bed: bed
            }
            .set { ch_phased_scatter_in }

        BCFTOOLS_VIEW_PHASING(
            ch_phased_scatter_in.vcf,
            ch_phased_scatter_in.bed,
            [],
            [],
        )
        ch_snv_vcf_for_annotation = BCFTOOLS_VIEW_PHASING.out.vcf
        ch_snv_index_for_annotation = BCFTOOLS_VIEW_PHASING.out.tbi
        ch_sv_vcf_for_annotation = PHASING.out.phased_family_svs
        ch_sv_index_for_annotation = PHASING.out.phased_family_svs_tbi
    }
    else {
        // Guarding against Nexflow trying to bind uninitialized channels even though we don't run annotation without SNVs
        if (!val_skip_snv_calling) {
            ch_snv_vcf_for_annotation = family_snv_vcf
            ch_snv_index_for_annotation = family_snv_index
        }
        if (!val_skip_sv_calling) {
            ch_sv_vcf_for_annotation = CALL_SVS.out.family_vcf
            ch_sv_index_for_annotation = CALL_SVS.out.family_tbi
        }
    }

    //
    // Annotate SNVs
    if (!val_skip_snv_annotation) {

        // Annotates family VCFs per variant call region
        ANNOTATE_SNVS(
            ch_snv_vcf_for_annotation,
            ch_echtvar_databases.map { _meta, databases -> databases }.collect(),
            ch_fasta,
            ch_fai,
            PREPARE_REFERENCES.out.vep_resources,
            val_vep_cache_version,
            ch_vep_plugin_files.collect(),
            ch_cadd_resources && ch_cadd_prescored_indels,
            val_echtvar_snv_databases,
            ch_cadd_header,
            ch_cadd_resources,
            ch_cadd_prescored_indels,
            val_pre_vep_snv_filter_expression != '',
        )

        ANNOTATE_SNVS.out.vcf
            .multiMap { meta, vcf ->
                clinical: [meta + [set: "clinical"], vcf]
                research: [meta + [set: "research"], vcf]
            }
            .set { ch_clin_research_snvs_vcf }
        ch_clin_research_snvs_vcf.clinical

        ch_clin_research_snvs_vcf.research.set { ch_ann_csq_pli_snv_in }

        if (val_filter_variants_hgnc_ids || val_filter_snvs_expression != '') {

            FILTER_VARIANTS_SNVS(
                ch_clin_research_snvs_vcf.clinical,
                ch_hgnc_ids,
                val_filter_snvs_expression,
                val_filter_variants_hgnc_ids,
            )

            ch_ann_csq_pli_snv_in = ch_ann_csq_pli_snv_in.mix(FILTER_VARIANTS_SNVS.out.vcf)
        }

        // This is really only required for ranking, could consider moving it there?
        ANN_CSQ_PLI_SNV(
            ch_ann_csq_pli_snv_in,
            ch_variant_consequences_snvs,
        )

        ANN_CSQ_PLI_SNV.out.vcf
            .join(ANN_CSQ_PLI_SNV.out.tbi, failOnMismatch: true, failOnDuplicate: true)
            .set { ch_vcf_tbi_per_region }
    }

    //
    // Concatenate and sort annotated SNVs for chromograph - requires an AF-tag, e.g. gnomad_af
    //
    def split_family_vcf_for_chromograph = !val_skip_chromograph && val_plot_chromograph_autozygosity && !val_skip_snv_annotation

    if (split_family_vcf_for_chromograph) {

        ANNOTATE_SNVS.out.vcf
            .join(ANNOTATE_SNVS.out.tbi, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, vcf, tbi ->
                def new_meta = [id: meta.family_id, num_intervals: meta.num_intervals]
                [groupKey(new_meta, new_meta.num_intervals), vcf, tbi]
            }
            .groupTuple()
            .set { ch_concat_sort_annotated_snvs_input }

        CONCAT_SORT_ANNOTATED_SNVS(
            ch_concat_sort_annotated_snvs_input
        )

        // Transpose family-level VCFs and add sample IDs by combining with samplesheet meta
        ch_samplesheet
            .map { meta, _files -> [id: meta.id, family_id: meta.family_id] }
            .unique()
            .combine(
                CONCAT_SORT_ANNOTATED_SNVS.out.vcf.join(CONCAT_SORT_ANNOTATED_SNVS.out.index, failOnMismatch: true, failOnDuplicate: true)
            )
            .filter { sample_info, vcf_meta, _vcf, _tbi ->
                sample_info.family_id == vcf_meta.id
            }
            .map { sample_info, _vcf_meta, vcf, tbi ->
                [sample_info, vcf, tbi]
            }
            .set { ch_bcftools_view_chromograph_input }

        BCFTOOLS_VIEW_CHROMOGRAPH(
            ch_bcftools_view_chromograph_input,
            [],
            [],
            [],
        )
    }

    if (!val_skip_chromograph) {
        CHROMOGRAPH(
            ch_bam_bai,
            split_family_vcf_for_chromograph ? BCFTOOLS_VIEW_CHROMOGRAPH.out.vcf : channel.empty(),
            split_family_vcf_for_chromograph ? BCFTOOLS_VIEW_CHROMOGRAPH.out.tbi : channel.empty(),
            val_plot_chromograph_coverage,
            val_plot_chromograph_autozygosity,
        )
    }

    //
    // Ranks family VCFs per variant call region
    // Can only run if samplesheet has affected samples
    //
    if (!val_skip_rank_variants) {

        // Create PED with updated sex - per family
        SOMALIER_PED_FAMILY(
            ch_bam.map { meta, _files -> [[id: meta.family_id], meta] }.groupTuple()
        )

        // Give PED file SNV meta so they can be joined later in the subworkflow.
        // Since we don't always have matching number of ped files and call regions
        // we need to combine and filter instead of join
        ANN_CSQ_PLI_SNV.out.vcf
            .map { meta, _vcf -> [[id: meta.family_id], meta] }
            .combine(SOMALIER_PED_FAMILY.out.ped)
            .filter { family_id_snv, _meta, family_id_ped, _ped -> family_id_snv == family_id_ped }
            .map { _family_id_snv, meta, _family_id_ped, ped -> [meta, ped] }
            .set { ch_snv_ranking_ped_file }

        // Only run if we have affected individuals
        RANK_VARIANTS_SNV(
            addChildWithTwoParentsToMeta(ANN_CSQ_PLI_SNV.out.vcf, ch_samplesheet, 'family_id'),
            addChildWithTwoParentsToMeta(ch_snv_ranking_ped_file, ch_samplesheet, 'family_id'),
            ch_genmod_reduced_penetrance,
            ch_genmod_score_config_snvs,
        )

        RANK_VARIANTS_SNV.out.vcf
            .join(RANK_VARIANTS_SNV.out.tbi, failOnMismatch: true, failOnDuplicate: true)
            .set { ch_vcf_tbi_per_region }
    }

    //
    // Concatenate and sort ranked SNVs, sort and publish
    //
    if (!val_skip_snv_calling) {
        ch_vcf_tbi_per_region
            .map { meta, vcf, tbi ->
                def new_meta = [id: meta.family_id, set: meta.set, sample_ids: meta.sample_ids, num_intervals: meta.num_intervals]
                [groupKey(new_meta, meta.num_intervals), vcf, tbi]
            }
            .groupTuple()
            .set { ch_concat_sort_input }

        CONCAT_SORT_RANKED_SNVS(
            ch_concat_sort_input
        )
    }

    //
    // Run Peddy
    //
    if (!val_skip_snv_calling && !val_skip_peddy) {

        CONCAT_SORT_RANKED_SNVS.out.vcf
            .join(CONCAT_SORT_RANKED_SNVS.out.index, failOnMismatch: true, failOnDuplicate: true)
            .filter { meta, _vcf, _tbi -> meta.set == "research" }
            .set { ch_peddy_in }

        PEDDY(
            ch_peddy_in,
            ch_samplesheet_pedfile,
            ch_peddy_sites,
        )
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped.map { _meta, metrics -> metrics }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.het_check_csv.map { _meta, metrics -> metrics }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped_check_csv.map { _meta, metrics -> metrics }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped_check_rel_difference_csv.map { _meta, metrics -> metrics }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.het_check_png.map { _meta, metrics -> metrics }.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped_check_png.map { _meta, metrics -> metrics }.collect().ifEmpty([]))
    }

    //
    // Annotate SVs
    //
    if (!val_skip_sv_annotation) {

        ANNOTATE_SVS(
            ch_sv_vcf_for_annotation,
            ch_fasta,
            ch_svdb_sv_databases,
            PREPARE_REFERENCES.out.vep_resources,
            val_vep_cache_version,
            ch_vep_plugin_files.collect(),
        )

        ANNOTATE_SVS.out.vcf
            .multiMap { meta, vcf ->
                clinical: [meta + [set: "clinical"], vcf]
                research: [meta + [set: "research"], vcf]
            }
            .set { ch_clin_research_svs_vcf }

        ch_clin_research_svs_vcf.research.set { ch_ann_csq_svs_in }

        //
        // Filter SVs
        //
        if (val_filter_variants_hgnc_ids || val_filter_svs_expression != '') {

            FILTER_VARIANTS_SVS(
                ch_clin_research_svs_vcf.clinical,
                ch_hgnc_ids,
                val_filter_svs_expression,
                val_filter_variants_hgnc_ids,
            )

            ch_ann_csq_svs_in = ch_ann_csq_svs_in.mix(FILTER_VARIANTS_SVS.out.vcf)
        }

        ANN_CSQ_PLI_SVS(
            ch_ann_csq_svs_in,
            ch_variant_consequences_svs,
        )
    }

    //
    // Rank SVs
    //
    if (!val_skip_rank_variants) {

        // Give PED file SVs meta so they can be joined later in the subworkflow.
        ANN_CSQ_PLI_SVS.out.vcf
            .combine(SOMALIER_PED_FAMILY.out.ped)
            .filter { vcf_meta, _vcf, ped_meta, _ped -> vcf_meta.id == ped_meta.id }
            .map { vcf_meta, _vcf, _ped_meta, ped -> [vcf_meta, ped] }
            .set { ch_sv_ranking_ped_file }

        RANK_VARIANTS_SVS(
            addChildWithTwoParentsToMeta(ANN_CSQ_PLI_SVS.out.vcf, ch_samplesheet, 'id'),
            addChildWithTwoParentsToMeta(ch_sv_ranking_ped_file, ch_samplesheet, 'id'),
            ch_genmod_reduced_penetrance,
            ch_genmod_score_config_svs,
        )
    }

    //
    // Collect and publish SVs
    //
    if (!val_skip_sv_calling) {

        ch_collect_svs = val_skip_sv_annotation
            ? ch_sv_vcf_for_annotation
            : val_skip_rank_variants
                ? ANN_CSQ_PLI_SVS.out.vcf
                : RANK_VARIANTS_SVS.out.vcf

        BCFTOOLS_VIEW_SV(
            ch_collect_svs.map { meta, vcf -> [meta, vcf, []] },
            [],
            [],
            [],
        )
    }

    //
    // Create methylation pileups with modkit or pbcpgtools, create methylation profile with methbat for pacbio
    //
    if (!val_skip_methylation_calling && val_run_modkit) {
        CALL_METHYLATION_MODKIT(
            !val_skip_phasing ? PHASING.out.haplotagged_bam_bai : ch_bam_bai,
            ch_fasta,
            ch_fai,
            ch_modkit_call_regions,
            val_bigwig_modcodes,
        )
    }

    if (!val_skip_methylation_calling && val_run_methbat) {
        CALL_METHYLATION_METHBAT(
            !val_skip_phasing ? PHASING.out.haplotagged_bam_bai : ch_bam_bai,
            ch_methbat_regions,
        )
    }

    //
    // Call repeat expansions with TRGT
    //
    if (!val_skip_repeat_calling) {
        if (val_str_caller == "trgt") {
            CALL_REPEAT_EXPANSIONS_TRGT(
                PHASING.out.haplotagged_bam_bai,
                ch_fasta,
                ch_fai,
                ch_str_bed,
                cram_output,
                ch_vcfexpress_prelude,
            )

            ch_repeat_expansions = CALL_REPEAT_EXPANSIONS_TRGT.out.family_vcf
        }
        else if (val_str_caller == "strdust") {
            CALL_REPEAT_EXPANSIONS_STRDUST(
                PHASING.out.haplotagged_bam_bai,
                ch_fasta,
                ch_fai,
                ch_str_bed,
                ch_vcfexpress_prelude,
            )
        }
    }

    //
    // Annotate repeat expansions with Stranger
    //
    if (!val_skip_repeat_annotation) {
        ANNOTATE_REPEAT_EXPANSIONS(
            ch_repeat_expansions,
            ch_strdrop_training_set_json,
            [[], []],
            ch_stranger_repeat_catalog,
            val_strdrop_training_set_json,
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${val_outdir}/pipeline_info",
            name: 'nallo_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )

    ch_multiqc_custom_methods_description = val_multiqc_methods_description
        ? file(val_multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )
    ch_methods_description_citation = citationBibliographyText(
        topic_versions_string,
        file("${projectDir}/assets/software_references.yml"),
        'citation',
    )
    ch_methods_description_bibliography = citationBibliographyText(
        topic_versions_string,
        file("${projectDir}/assets/software_references.yml"),
        'bibliography',
    )
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    // sort: false // preserve order for correct yaml structure
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.concat(ch_methods_description_citation).concat(ch_methods_description_bibliography).flatten().collectFile(
            name: 'methods_description_mqc.yaml',
            sort: false,
        )
    )

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'multiqc'],
                files,
                val_multiqc_config
                    ? file(val_multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                val_multiqc_logo ? file(val_multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    aligned_assemblies = val_skip_genome_assembly ? channel.empty() : cram_output ? ALIGN_ASSEMBLIES.out.cram.join(ALIGN_ASSEMBLIES.out.crai) : ALIGN_ASSEMBLIES.out.bam.join(ALIGN_ASSEMBLIES.out.bai) // channel: [ val(meta), path(bam/cram), path(bai/crai) ]

}

// Check if a family has a child with two parents,
// and add this information to the input variant channel meta as 'child_with_two_parents_in_family'.
// This is used to determine compound ranking thresholds and penalties in genmod.
def addChildWithTwoParentsToMeta(input, samplesheet, family_id_key) {
    samplesheet
        .map { meta, _files ->
            [meta.family_id, meta]
        }
        .groupTuple()
        .map { family_id, metas ->
            [id: family_id, child_with_two_parents_in_family: metas.any { meta -> meta.two_parents }]
        }
        .combine(input)
        .filter { family_meta, vcf_meta, _file -> vcf_meta[family_id_key] == family_meta.id }
        .map { family_meta, vcf_meta, file ->
            def new_meta = vcf_meta + [child_with_two_parents_in_family: family_meta.child_with_two_parents_in_family]
            [new_meta, file]
        }
}
