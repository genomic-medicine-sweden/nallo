include { samplesheetToList                                      } from 'plugin/nf-schema'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGN                                                  } from '../subworkflows/local/align'
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
include { CALL_MITOCHONDRIAL_VARIANTS                            } from '../subworkflows/local/call_mitochondrial_variants'
include { CALL_PARALOGS                                          } from '../subworkflows/local/call_paralogs'
include { CALL_REPEAT_EXPANSIONS_STRDUST                         } from '../subworkflows/local/call_repeat_expansions_strdust'
include { CALL_REPEAT_EXPANSIONS_TRGT                            } from '../subworkflows/local/call_repeat_expansions_trgt'
include { CALL_SNVS                                              } from '../subworkflows/local/call_snvs'
include { CALL_SVS                                               } from '../subworkflows/local/call_svs'
include { GENOME_ASSEMBLY                                        } from '../subworkflows/local/genome_assembly'
include { GVCF_GLNEXUS_NORM_VARIANTS                             } from '../subworkflows/local/gvcf_glnexus_norm_variants'
include { CALL_METHYLATION_MODKIT                                } from '../subworkflows/local/call_methylation_modkit'
include { CALL_METHYLATION_METHBAT                               } from '../subworkflows/local/call_methylation_methbat'
include { MERGE_SVS                                              } from '../subworkflows/local/merge_svs/main'
include { PHASING                                                } from '../subworkflows/local/phasing'
include { PREPARE_GENS_INPUTS                                    } from '../subworkflows/local/prepare_gens_inputs'
include { PREPARE_REFERENCES                                     } from '../subworkflows/local/prepare_references'
include { QC_ALIGNED_READS                                       } from '../subworkflows/local/qc_aligned_reads'
include { QC_SNVS                                                } from '../subworkflows/local/qc_snvs'
include { RANK_VARIANTS                                          } from '../subworkflows/local/rank_variants'
include { SCATTER_GENOME                                         } from '../subworkflows/local/scatter_genome'
include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP as FILTER_VARIANTS_SNVS } from '../subworkflows/nf-core/vcf_filter_bcftools_ensemblvep/main'
include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP as FILTER_VARIANTS_SVS  } from '../subworkflows/nf-core/vcf_filter_bcftools_ensemblvep/main'
include { VCF_CONCAT_NORM_VARIANTS                               } from '../subworkflows/local/vcf_concat_norm_variants'
include { VCF_CONCAT_SORT_VARIANTS as CONCAT_SORT_ANNOTATED_SNVS } from '../subworkflows/local/vcf_concat_sort_variants/main'
include { VCF_CONCAT_SORT_VARIANTS as CONCAT_SORT_RANKED_SNVS    } from '../subworkflows/local/vcf_concat_sort_variants/main'
include { VCF_CONCAT_SORT_VARIANTS as CONCAT_SORT_GENS           } from '../subworkflows/local/vcf_concat_sort_variants/main'
include { VCF_CONCAT_SORT_VARIANTS as CONCAT_SORT_PEDDY          } from '../subworkflows/local/vcf_concat_sort_variants/main'
include { ANNOTATE_METHYLATION                                   } from '../subworkflows/local/annotate_methylation'
include { PORTELLO                                               } from '../subworkflows/local/portello/main'

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
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_MITO_SNVS           } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CHROMOGRAPH             } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SV                      } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_PHASING                 } from '../modules/nf-core/bcftools/view/main'
include { MINIMAP2_ALIGN                                         } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                                         } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                                         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_CONVERT                                       } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CALMD                                         } from '../modules/nf-core/samtools/calmd/main'
include { MULTIQC                                                } from '../modules/nf-core/multiqc/main'
include { PEDDY                                                  } from '../modules/nf-core/peddy/main'
include { SPLITUBAM                                              } from '../modules/nf-core/splitubam/main'
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
    ch_fai
    ch_genmod_reduced_penetrance
    ch_genmod_score_config_snvs
    ch_genmod_score_config_svs
    ch_gens_baf_positions
    ch_gens_coverage_bins
    ch_gens_panel_of_normals_female
    ch_gens_panel_of_normals_male
    ch_glnexus_config
    ch_hgnc_ids
    ch_samplesheet
    ch_cramino_regions
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
    val_alignment_processes
    val_assembly_aligner
    val_bigwig_modcodes
    val_convert_unphased_aligned_reads_to_cram
    val_cram_output
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
    val_mosdepth_regions
    val_cramino_regions
    val_multiqc_config
    val_multiqc_logo
    val_multiqc_methods_description
    val_outdir
    val_paraphrase_output_format
    val_phaser
    val_whatshap_pedigree_phasing
    val_plot_chromograph_autozygosity
    val_plot_chromograph_coverage
    val_pre_vep_snv_filter_expression
    val_premapped
    val_read_aligner
    val_sentieon_tech
    val_skip_alignment
    val_skip_annotate_paralogs
    val_skip_call_paralogs
    val_skip_chromograph
    val_skip_genome_assembly
    val_skip_methbat
    val_skip_methylation_annotation
    val_skip_modkit
    val_skip_peddy
    val_skip_phasing
    val_skip_portello
    val_skip_prepare_gens_input
    val_skip_qc
    val_skip_rank_variants
    val_skip_repeat_annotation
    val_skip_mitochondrial_calling
    val_skip_sambamba_depth
    val_skip_sex_check
    val_skip_snv_annotation
    val_skip_snv_calling
    val_skip_trgt
    val_skip_sv_annotation
    val_skip_sv_calling
    val_skip_strdust
    val_snv_caller
    val_snv_calling_processes
    val_snv_call_regions
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
            ch_fai,
            ch_vep_cache_unprocessed,
            val_fasta.endsWith('.gz'),
            val_vep_cache && val_vep_cache.endsWith("tar.gz"),
        )

        // Gather indices
        ch_fasta = PREPARE_REFERENCES.out.fasta
        ch_fai = PREPARE_REFERENCES.out.fai
    }


    /*
     * Map reads to reference
     */
    if (!val_skip_alignment) {

        if (!val_premapped) {

            CONVERT_INPUT_FASTQS(
                ch_samplesheet,
                false,
                true,
            )

            if (val_alignment_processes > 1) {
                SPLITUBAM(CONVERT_INPUT_FASTQS.out.bam)
                ch_unmapped = SPLITUBAM.out.bam.transpose()
            }
            else {
                ch_unmapped = CONVERT_INPUT_FASTQS.out.bam
            }
        }
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
        // Since starting with FASTQs is a rare case, no splitting of FASTQs alone just for the assembly is implemented

        CONVERT_INPUT_BAMS(
            val_skip_alignment || val_premapped || (val_alignment_processes == 1) ? ch_samplesheet : SPLITUBAM.out.bam.transpose(),
            true,
            false,
        )

        // contains all FASTQ files, including those not converted
        ch_genome_assembly_input = CONVERT_INPUT_BAMS.out.fastq.groupTuple()

        // Hifiasm assembly
        // Concatenate haplotypes per sample if portello is not skipped so it can be used for reads-to-assembly alignment
        GENOME_ASSEMBLY(
            ch_genome_assembly_input,
            val_hifiasm_mode == "trio-binning",
            !val_skip_portello,
        )

        ALIGN_ASSEMBLIES(
            GENOME_ASSEMBLY.out.assembled_haplotypes,
            ch_fasta,
            ch_fai,
            val_cram_output,
            val_assembly_aligner,
        )

        ch_assembly_bam_bai = ALIGN_ASSEMBLIES.out.unfiltered_bam.join(ALIGN_ASSEMBLIES.out.unfiltered_bai, failOnMismatch: true, failOnDuplicate: true)
    }

    if (!val_skip_alignment) {

        if (!val_premapped) {
            /*
             * Create a grouping key per sample that records the number of split files,
             * allowing downstream merging to trigger as soon as all alignments of a sample are ready.
             */
            ch_reads_grouping_key = ch_unmapped
                .groupTuple()
                .map { meta, files -> tuple(meta.id, files.size()) }

            // Add original file name to meta to join correct alignments and indexes
            ch_align_in = ch_unmapped.map { meta, bam -> tuple(meta + [file: bam.name], bam) }

            // If portello is not skipped, we need to align the reads to the concatenated haplotypes from the assembly,
            // otherwise we can align to the reference genome.
            ALIGN(
                ch_align_in,
                !val_skip_portello ? GENOME_ASSEMBLY.out.concatenated_haplotypes : ch_fasta,
                val_read_aligner,
                val_skip_portello,
            )

            ch_aligned_for_merge = ALIGN.out.bam
                .join(ALIGN.out.index, failOnMismatch: true, failOnDuplicate: true)
                .combine(ch_reads_grouping_key)
                .filter { bam_meta, _bam, _bai, group_id, _group_size ->
                    bam_meta.id == group_id
                }
                .map { bam_meta, bam, bai, _group_id, group_size ->
                    tuple(groupKey(bam_meta - bam_meta.subMap('file'), group_size), bam, bai)
                }
                .groupTuple()
                .map { key, bams, bais -> tuple(key.getGroupTarget(), bams, bais) }
                .map { meta, bams, bais ->
                    // Keep BAM and BAI pairing while enforcing deterministic order.
                    def bam_bai_pairs = [bams, bais]
                        .transpose()
                        .sort { left, right ->
                            left[0].getName() <=> right[0].getName()
                        }
                    [meta, bam_bai_pairs.collect { pair -> pair[0] }, bam_bai_pairs.collect { pair -> pair[1] }]
                }
        }
        else {

            // If bams are premapped, just merge them (ONT machines output several BAMs per sample)
            // SAMTOOLS_MERGE expects indexes in the input but is happy to merge them if the indexes are missing
            ch_aligned_for_merge = ch_samplesheet
                .groupTuple()
                .map { meta, reads -> [meta, reads, []] }
        }

        SAMTOOLS_MERGE(
            ch_aligned_for_merge,
            [[], [], [], []],
        )

        ch_aligned_bam = SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_MERGE.out.index, failOnMismatch: true, failOnDuplicate: true)

        // Publish alignments as CRAM if requested
        if (val_convert_unphased_aligned_reads_to_cram) {
            SAMTOOLS_CONVERT(
                ch_aligned_bam,
                ch_fasta.join(ch_fai).collect(),
            )
        }

        if (!val_skip_portello) {

            PORTELLO(
                ch_aligned_bam,
                ch_assembly_bam_bai,
                ch_fasta,
            )

            ch_aligned_bam = PORTELLO.out.bam.join(PORTELLO.out.bai, failOnMismatch: true, failOnDuplicate: true)
        }

        if (val_sv_callers_to_run.contains("sniffles") && val_read_aligner == "pbmm2") {
            SAMTOOLS_CALMD(
                val_skip_portello ? SAMTOOLS_MERGE.out.bam : PORTELLO.out.bam,
                ch_fasta.join(ch_fai).collect(),
            )

            SAMTOOLS_INDEX(SAMTOOLS_CALMD.out.bam)

            ch_aligned_bam = SAMTOOLS_CALMD.out.bam.join(SAMTOOLS_INDEX.out.index, failOnMismatch: true, failOnDuplicate: true)
        }

        //
        // Create PED from samplesheet
        //
        ch_samplesheet_ped_in = ch_samplesheet
            .map { meta, _files -> [[id: meta.project], meta] }
            .groupTuple()

        SAMPLESHEET_PED(ch_samplesheet_ped_in)

        ch_samplesheet_pedfile = SAMPLESHEET_PED.out.ped.collect()

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

        if (!val_skip_rank_variants || !val_skip_phasing) {
            // Create PED files with updated (infered sex) per family
            SOMALIER_PED_FAMILY(
                ch_bam.map { meta, _files -> [[id: meta.family_id], meta] }.groupTuple()
            )
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
            val_cramino_regions,
            ch_cramino_regions,
        )
        ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNED_READS.out.fastqc_zip.collect { _meta, metrics -> metrics }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNED_READS.out.mosdepth_summary.collect { _meta, metrics -> metrics })
        ch_multiqc_files = ch_multiqc_files.mix(
            (val_mosdepth_regions
                ? QC_ALIGNED_READS.out.mosdepth_regions_dist
                : QC_ALIGNED_READS.out.mosdepth_global_dist).collect { _meta, metrics -> metrics }
        )
    }

    /*
     * Call paralogous genes with paraphase
     */
    if (!val_skip_call_paralogs) {
        CALL_PARALOGS(
            ch_bam_bai,
            ch_fasta,
            ch_fai,
            val_cram_output,
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

        ch_bed_intervals = SCATTER_GENOME.out.bed_nuclear_intervals.map { meta, bed, num_intervals -> [meta + [caller: val_snv_caller], bed, num_intervals] }
        ch_mitochondrial_bed = SCATTER_GENOME.out.bed_mitochondrial_intervals.map { meta, bed, _num_intervals -> [meta, bed] }

        def ch_num_intervals = ch_bed_intervals.map { _meta, _bed, num_intervals -> num_intervals }.first()

        if (!val_skip_mitochondrial_calling) {
            CALL_MITOCHONDRIAL_VARIANTS(
                ch_bam_bai,
                ch_fasta,
                ch_fai,
                ch_par,
                ch_mitochondrial_bed,
                val_mitochondrial_caller,
            )

            /*
             * The meta.caller is needed for GVCF_GLNEXUS_NORM_VARIANTS workflow to process the VCFs differently.
             * the number of intervals is used in groupKey downstream, num_intervals should be the same in the nuclear and mitochondrial channels.
             */
            ch_mitochondrial = CALL_MITOCHONDRIAL_VARIANTS.out.mitochondrial_snv_vcf
                .join(CALL_MITOCHONDRIAL_VARIANTS.out.mitochondrial_snv_tbi, failOnMismatch: true, failOnDuplicate: true)
                .combine(ch_num_intervals)
                .multiMap { meta, vcf, tbi, num_intervals ->
                    vcf: [meta + [caller: val_mitochondrial_caller, genome: 'mitochondrial', num_intervals: num_intervals + 1], vcf]
                    index: [meta + [caller: val_mitochondrial_caller, genome: 'mitochondrial', num_intervals: num_intervals + 1], tbi]
                }
        }
        else {
            ch_mitochondrial = channel.empty()
                .multiMap { it ->
                    vcf: it
                    index: it
                }
        }

        // Combine the BED intervals with BAM/BAI files to create a region-bam-bai for each sample.
        // This uses the whole BAM files for each region instead of splitting them.
        ch_call_snvs_input = ch_bam_bai
            .combine(ch_bed_intervals)
            .map { meta, bam, bai, bed_meta, bed, num_intervals ->
                [meta + [genome: bed_meta.genome, num_intervals: num_intervals, region: bed], bam, bai, bed]
            }

        CALL_SNVS(
            ch_call_snvs_input,
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

        /*
         * Group and create (GVCF_GLNEXUS_NORM_VARIANTS) a merged and normalized VCF, containing one region with all samples, to be used in annotation and ranking.
         * We add the total number of intervals for grouping later (account for mitochondrial interval).
         * Done for CALL_SNVS.out.gvcf but CALL_SNVS.out.vcf is used in QC_SNVS subworkflow, for now we do not want the mitochondrial variants to be included
         * as it will make the deepvariant report harder to interpret.
         */
        def variants_to_merge_per_family = CALL_SNVS.out.gvcf
            .join(CALL_SNVS.out.gvcf_index, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, gvcf, index ->
                def num_intervals = val_skip_mitochondrial_calling ? meta.num_intervals : meta.num_intervals + 1
                [[id: meta.region.name, family_id: meta.family_id, genome: meta.genome, num_intervals: num_intervals, caller: val_snv_caller], gvcf, index]
            }
            .mix(
                ch_mitochondrial.vcf.join(ch_mitochondrial.index, failOnMismatch: true, failOnDuplicate: true).map { meta, vcf, tbi ->
                    [[id: meta.genome, family_id: meta.family_id, genome: meta.genome, num_intervals: meta.num_intervals, caller: meta.caller], vcf, tbi]
                }
            )
            .groupTuple()
            .multiMap { meta, gvcfs, indexes ->
                gvcf: [meta, gvcfs]
                index: [meta, indexes]
            }

        // SCATTER_GENOME.out.bed contains all regions, but we could probably pass the region BED that actually matches the variants instead...
        GVCF_GLNEXUS_NORM_VARIANTS(
            variants_to_merge_per_family.gvcf,
            variants_to_merge_per_family.index,
            SCATTER_GENOME.out.bed,
            ch_fasta,
            ch_fai,
            ch_vcfexpress_prelude,
            ch_glnexus_config,
        )

        // Grouping VCF, containing one sample with all regions except chrM, as we do not want mitochondrial variants in the deepvariant report for now.
        ch_variants_to_concat_per_sample = CALL_SNVS.out.vcf
            .map { meta, vcf ->
                def new_meta = meta - meta.subMap('region', 'genome')
                [groupKey(new_meta, new_meta.num_intervals), vcf]
            }
            .groupTuple()
            .map { meta, vcfs ->
                [meta - meta.subMap('num_intervals'), vcfs]
            }

        // Create a concatenated and normalized VCF, containing one sample with all regions.
        VCF_CONCAT_NORM_VARIANTS(
            ch_variants_to_concat_per_sample,
            ch_fasta,
            val_snv_caller,
            ch_vcfexpress_prelude,
        )

        // These contains RefCalls
        sample_snv_vcf = VCF_CONCAT_NORM_VARIANTS.out.vcf
        sample_snv_index = VCF_CONCAT_NORM_VARIANTS.out.index


        // SNV QC
        // Can we use the normalized VCF here, for DV vcfstatsreport?
        QC_SNVS(
            VCF_CONCAT_NORM_VARIANTS.out.bcftools_concat_vcf.map { meta, vcf -> [meta - meta.subMap('caller'), vcf] },
            sample_snv_vcf,
            sample_snv_index,
            val_snv_caller.equals("deepvariant"),
        )
        ch_multiqc_files = ch_multiqc_files.mix(QC_SNVS.out.stats.collect { _meta, metrics -> metrics }.ifEmpty([]))

        // Set family_snv_vcf and family_snv_index for clarity
        family_snv_vcf = GVCF_GLNEXUS_NORM_VARIANTS.out.vcf
        family_snv_index = GVCF_GLNEXUS_NORM_VARIANTS.out.index

        ch_snvs_per_family_unannotated_vcf_tbi = family_snv_vcf.join(family_snv_index, failOnMismatch: true, failOnDuplicate: true)
    }

    if (!val_skip_prepare_gens_input) {
        ch_gvcfs = CALL_SNVS.out.gvcf
            .join(CALL_SNVS.out.gvcf_index)
            .map { meta, gvcf, gvcf_index ->
                def sample_meta = meta - meta.subMap(['region', 'num_intervals', 'genome'])
                [sample_meta, gvcf, gvcf_index]
            }
            .groupTuple()

        CONCAT_SORT_GENS(
            ch_gvcfs
        )

        ch_gvcf_tbi = CONCAT_SORT_GENS.out.vcf.join(CONCAT_SORT_GENS.out.index)

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
            ch_sv_call_regions,
            val_sv_call_regions,
            val_force_sawfish_joint_call_single_samples,
            val_create_hificnv_maf_track,
            val_create_sawfish_maf_track,
        )

        MERGE_SVS(
            CALL_SVS.out.vcf,
            val_sv_callers_to_merge.split(',').collect { caller -> caller.toLowerCase().trim() },
            val_sv_callers_merge_priority.split(',').collect { caller -> caller.toLowerCase().trim() },
            ch_vcfexpress_prelude,
        )
    }

    //
    // Phase SNVs, SVs and INDELs
    //
    if (!val_skip_phasing) {

        ch_family_to_samples = ch_samplesheet
            .map { meta, _files -> [meta.family_id, meta.id] }
            .groupTuple()
            .map { family_id, sample_ids ->
                [[id: family_id], sample_ids.unique()]
            }

        /*
         * The VCFs are split by calling regions but we need whole-genome VCFs for phasing, we first group by family and then concatenate the VCFs of the same family together.
         * We group only nuclear vcf for phasing as phasing mitochondrial variants is not relevant."num_intervals - 1" happens because groupKey should not expect the mitochondrial interval.
         */
        ch_bcftools_concat_phasing_in = family_snv_vcf
            .join(family_snv_index, failOnMismatch: true, failOnDuplicate: true)
            .filter { meta, _vcf, _tbi -> meta.genome == "nuclear" }
            .map { meta, vcf, tbi ->
                def new_meta = [id: meta.family_id, num_intervals: val_skip_mitochondrial_calling ? meta.num_intervals : meta.num_intervals - 1]
                [groupKey(new_meta, new_meta.num_intervals), vcf, tbi]
            }
            .groupTuple()
            .map { key, vcfs, tbis ->
                [key.getGroupTarget(), vcfs, tbis]
            }
            .map { meta, vcfs, tbis ->
                [meta - meta.subMap('num_intervals'), vcfs, tbis]
            }

        BCFTOOLS_CONCAT_PHASING(
            ch_bcftools_concat_phasing_in
        )

        // Provide a PED file to let whatshap activate pedigree phasing
        // Or pass 'empty_PED' if 'whatshap_pedigree_phasing==false'
        ch_ped_family = SOMALIER_PED_FAMILY.out.ped.map { meta, ped -> [[id: meta.id], val_whatshap_pedigree_phasing ? ped : []] }

        // Input is one VCF per family with all the regions (except chrM) and all the samples in the family
        PHASING(
            BCFTOOLS_CONCAT_PHASING.out.vcf,
            BCFTOOLS_CONCAT_PHASING.out.tbi,
            val_skip_sv_calling ? channel.empty() : MERGE_SVS.out.family_vcf,
            val_skip_sv_calling ? channel.empty() : MERGE_SVS.out.family_tbi,
            ch_bam_bai,
            ch_family_to_samples,
            ch_fasta,
            ch_fai,
            val_phaser,
            !val_skip_sv_calling,
            val_cram_output,
            ch_ped_family,
        )

        ch_multiqc_files = ch_multiqc_files.mix(PHASING.out.stats.collect { _meta, txt -> txt }.ifEmpty([]))

        // Scatter nuclear whole-genome phased SNV VCFs back into regions for annotation.
        ch_phased_scatter_in = PHASING.out.phased_family_snvs
            .join(PHASING.out.phased_family_snvs_tbi, failOnMismatch: true, failOnDuplicate: true)
            .combine(ch_bed_intervals)
            .multiMap { vcf_meta, vcf, tbi, bed_meta, bed, num_intervals ->
                vcf: [vcf_meta + [id: bed.name, family_id: vcf_meta.id, genome: bed_meta.genome, num_intervals: num_intervals], vcf, tbi]
                bed: bed
            }

        // Split the vcf into regions for annotation
        BCFTOOLS_VIEW_PHASING(
            ch_phased_scatter_in.vcf,
            ch_phased_scatter_in.bed,
            [],
            [],
        )
        // Add the mitochondrial VCF after phasing SNVs.
        // Add +1 to the num_intervals of nuclear channel to account for mitochondrial region
        ch_snv_vcf_tbi_nuclear_for_annotation = BCFTOOLS_VIEW_PHASING.out.vcf
            .join(BCFTOOLS_VIEW_PHASING.out.tbi, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, vcf, tbi -> [meta + [num_intervals: val_skip_mitochondrial_calling ? meta.num_intervals : meta.num_intervals + 1], vcf, tbi] }

        ch_snv_vcf_tbi_mitochondrial_for_annotation = ch_snvs_per_family_unannotated_vcf_tbi.filter { meta, _vcf, _tbi -> meta.genome == "mitochondrial" }

        ch_snv_vcf_tbi_nuclear_mitochondrial_for_annotation = ch_snv_vcf_tbi_nuclear_for_annotation
            .mix(ch_snv_vcf_tbi_mitochondrial_for_annotation)
            .multiMap { meta, vcf, tbi ->
                vcf: [meta, vcf]
                tbi: [meta, tbi]
            }

        // Set phased SVs for annotation if SV calling is not skipped
        ch_sv_vcf_for_annotation = PHASING.out.phased_family_svs
        ch_sv_index_for_annotation = PHASING.out.phased_family_svs_tbi
    }
    else {
        ch_snv_vcf_tbi_nuclear_mitochondrial_for_annotation = val_skip_snv_calling
            ? channel.empty()
            : family_snv_vcf.join(family_snv_index, failOnMismatch: true, failOnDuplicate: true).multiMap { meta, vcf, tbi ->
                vcf: [meta, vcf]
                index: [meta, tbi]
            }
        ch_sv_vcf_for_annotation = val_skip_sv_calling ? channel.empty() : MERGE_SVS.out.family_vcf
        ch_sv_index_for_annotation = val_skip_sv_calling ? channel.empty() : MERGE_SVS.out.family_tbi
    }

    // Annotate SNVs
    if (!val_skip_snv_annotation) {

        // Annotates family VCFs per variant call region
        ANNOTATE_SNVS(
            ch_snv_vcf_tbi_nuclear_mitochondrial_for_annotation.vcf,
            ch_echtvar_databases.map { _meta, databases -> databases }.collect().map { dbs -> [[id: 'echtvar_db'], dbs] },
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

        ch_clin_research_snvs_vcf = ANNOTATE_SNVS.out.vcf.multiMap { meta, vcf ->
            clinical: [meta + [set: "clinical"], vcf]
            research: [meta + [set: "research"], vcf]
        }

        ch_ann_csq_pli_snv_in = ch_clin_research_snvs_vcf.research

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

        ch_snvs_per_family_annotated_vcf_tbi = ANN_CSQ_PLI_SNV.out.vcf.join(ANN_CSQ_PLI_SNV.out.tbi, failOnMismatch: true, failOnDuplicate: true)
    }

    //
    // Concatenate and sort annotated SNVs for chromograph - requires an AF-tag, e.g. gnomad_af
    //
    def split_family_vcf_for_chromograph = !val_skip_chromograph && val_plot_chromograph_autozygosity && !val_skip_snv_annotation

    if (split_family_vcf_for_chromograph || (!val_skip_peddy && !val_skip_snv_annotation)) {

        ch_concat_sort_annotated_snvs_input = ANNOTATE_SNVS.out.unfiltered_vcf
            .join(ANNOTATE_SNVS.out.unfiltered_tbi, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, vcf, tbi ->
                def new_meta = [id: meta.family_id, num_intervals: meta.num_intervals]
                [groupKey(new_meta, new_meta.num_intervals), vcf, tbi]
            }
            .groupTuple()

        CONCAT_SORT_ANNOTATED_SNVS(
            ch_concat_sort_annotated_snvs_input
        )
    }

    if (split_family_vcf_for_chromograph) {
        // Transpose family-level VCFs and add sample IDs by combining with samplesheet meta
        ch_bcftools_view_chromograph_input = ch_samplesheet
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

        BCFTOOLS_VIEW_CHROMOGRAPH(
            ch_bcftools_view_chromograph_input,
            [],
            [],
            [],
        )
    }

    //
    // Run Peddy
    //
    if (!val_skip_peddy) {

        if (!val_skip_snv_annotation) {
            // Use already concatenated VCFs
            ch_peddy_in = CONCAT_SORT_ANNOTATED_SNVS.out.vcf.join(CONCAT_SORT_ANNOTATED_SNVS.out.index, failOnMismatch: true, failOnDuplicate: true)
        }
        else {
            // If we did not annotate, we did not concatenate the VCFs before, so we need to do that here.
            ch_concat_sort_peddy_in = ch_snvs_per_family_unannotated_vcf_tbi
                .map { meta, vcf, tbi -> [groupKey([id: meta.family_id], meta.num_intervals), vcf, tbi] }
                .groupTuple()
                .map { key, vcfs, tbis -> [key.getGroupTarget(), vcfs, tbis] }

            CONCAT_SORT_PEDDY(
                ch_concat_sort_peddy_in
            )

            ch_peddy_in = CONCAT_SORT_PEDDY.out.vcf.join(CONCAT_SORT_PEDDY.out.index, failOnMismatch: true, failOnDuplicate: true)
        }

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

        ch_clin_research_svs_vcf = ANNOTATE_SVS.out.vcf.multiMap { meta, vcf ->
            clinical: [meta + [set: "clinical"], vcf]
            research: [meta + [set: "research"], vcf]
        }

        ch_ann_csq_svs_in = ch_clin_research_svs_vcf.research

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

    /*
     * Ranks family VCFs per variant call region for SNVs,
     * and per family for SVs, based on provided ranking config.
     * Can only run if samplesheet has affected samples.
     */
    if (!val_skip_rank_variants) {

        ch_snvs_to_rank = buildRankVariantsInputChannel(
            ANN_CSQ_PLI_SNV.out.vcf,
            SOMALIER_PED_FAMILY.out.ped,
            'snv',
            ch_genmod_score_config_snvs,
            ch_samplesheet,
        )

        ch_svs_to_rank = buildRankVariantsInputChannel(
            ANN_CSQ_PLI_SVS.out.vcf.map { meta, vcf -> [meta + [family_id: meta.id], vcf] },
            SOMALIER_PED_FAMILY.out.ped,
            'sv',
            ch_genmod_score_config_svs,
            ch_samplesheet,
        )

        ch_rank_variants_input = ch_snvs_to_rank
            .mix(ch_svs_to_rank)
            .multiMap { meta, vcf, ped, score_config ->
                vcf: [meta, vcf]
                ped: [meta, ped]
                score_config: [meta, score_config]
            }

        RANK_VARIANTS(
            ch_rank_variants_input.vcf,
            ch_rank_variants_input.ped,
            ch_genmod_reduced_penetrance,
            ch_rank_variants_input.score_config,
        )

        ch_ranked_variants = RANK_VARIANTS.out.vcf
            .join(RANK_VARIANTS.out.tbi, failOnMismatch: true, failOnDuplicate: true)
            .branch { meta, _vcf, _tbi ->
                snvs: meta.variant_type == "snv"
                sv: meta.variant_type == "sv"
            }

        ch_snvs_per_family_ranked_vcf_tbi = ch_ranked_variants.snvs
    }

    //
    // Concatenate and sort SNVs, sort and publish
    //
    if (!val_skip_snv_calling) {
        def ch_snvs_per_family_to_concatenate = val_skip_rank_variants
            ? (val_skip_snv_annotation ? ch_snvs_per_family_unannotated_vcf_tbi : ch_snvs_per_family_annotated_vcf_tbi)
            : ch_snvs_per_family_ranked_vcf_tbi

        ch_concat_sort_input = ch_snvs_per_family_to_concatenate
            .map { meta, vcf, tbi ->
                def new_meta = [id: meta.family_id, set: meta.set, sample_ids: meta.sample_ids, num_intervals: meta.num_intervals]
                [groupKey(new_meta, meta.num_intervals), vcf, tbi]
            }
            .groupTuple()

        CONCAT_SORT_RANKED_SNVS(
            ch_concat_sort_input
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
                : ch_ranked_variants.sv.map { meta, vcf, _tbi -> [meta, vcf] }

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
    if (!val_skip_modkit) {
        CALL_METHYLATION_MODKIT(
            !val_skip_phasing ? PHASING.out.haplotagged_bam_bai : ch_bam_bai,
            ch_fasta,
            ch_fai,
            ch_modkit_call_regions,
            val_bigwig_modcodes,
        )
    }

    if (!val_skip_methbat) {
        CALL_METHYLATION_METHBAT(
            !val_skip_phasing ? PHASING.out.haplotagged_bam_bai : ch_bam_bai,
            ch_methbat_regions,
        )

        ch_methylation_profiles = CALL_METHYLATION_METHBAT.out.region_profile
    }

    if (!val_skip_methylation_annotation && !val_skip_methbat) {
        ANNOTATE_METHYLATION(
            ch_methylation_profiles
        )
    }

    //
    // Call repeat expansions with TRGT or strdust
    //
    if (!val_skip_trgt) {
        CALL_REPEAT_EXPANSIONS_TRGT(
            PHASING.out.haplotagged_bam_bai,
            ch_fasta,
            ch_fai,
            ch_str_bed,
            val_cram_output,
            ch_vcfexpress_prelude,
        )

        ch_repeat_expansions = CALL_REPEAT_EXPANSIONS_TRGT.out.family_vcf
    }
    if (!val_skip_strdust) {
        CALL_REPEAT_EXPANSIONS_STRDUST(
            PHASING.out.haplotagged_bam_bai,
            ch_fasta,
            ch_fai,
            ch_str_bed,
            ch_vcfexpress_prelude,
        )
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

    ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${val_outdir}/pipeline_info",
            name: 'nallo_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // MODULE: MultiQC
    //
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )

    def ch_multiqc_custom_methods_description = val_multiqc_methods_description
        ? file(val_multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )
    def ch_methods_description_citation = citationBibliographyText(
        topic_versions_string,
        file("${projectDir}/assets/software_references.yml"),
        'citation',
    )
    def ch_methods_description_bibliography = citationBibliographyText(
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
                [id: 'nallo'],
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
    aligned_assemblies_bai              = val_skip_genome_assembly ? channel.empty() : ALIGN_ASSEMBLIES.out.bai // channel: [ val(meta), path(bai) ]
    aligned_assemblies_bam              = val_skip_genome_assembly ? channel.empty() : ALIGN_ASSEMBLIES.out.bam // channel: [ val(meta), path(bam) ]
    aligned_assemblies_crai             = val_skip_genome_assembly ? channel.empty() : ALIGN_ASSEMBLIES.out.crai // channel: [ val(meta), path(crai) ]
    aligned_assemblies_cram             = val_skip_genome_assembly ? channel.empty() : ALIGN_ASSEMBLIES.out.cram // channel: [ val(meta), path(cram) ]
    aligned_assemblies_remapped_bam     = val_skip_portello ? channel.empty() : PORTELLO.out.bam // channel: [ val(meta), path(bam) ]
    aligned_assemblies_remapped_bai     = val_skip_portello ? channel.empty() : PORTELLO.out.bai // channel: [ val(meta), path(bai) ]
    aligned_haplotagged_reads_bai       = val_skip_phasing ? channel.empty() : PHASING.out.haplotagged_bam_bai.map { meta, _bam, bai -> [meta, bai] } // channel: [ val(meta), path(bai) ]
    aligned_haplotagged_reads_bam       = val_skip_phasing ? channel.empty() : PHASING.out.haplotagged_bam_bai.map { meta, bam, _bai -> [meta, bam] } // channel: [ val(meta), path(bam) ]
    aligned_haplotagged_reads_crai      = val_skip_phasing ? channel.empty() : PHASING.out.haplotagged_cram_crai.map { meta, _cram, crai -> [meta, crai] } // channel: [ val(meta), path(crai) ]
    aligned_haplotagged_reads_cram      = val_skip_phasing ? channel.empty() : PHASING.out.haplotagged_cram_crai.map { meta, cram, _crai -> [meta, cram] } // channel: [ val(meta), path(cram) ]
    aligned_reads_bai                   = val_skip_alignment ? channel.empty() : ch_aligned_bam.map { meta, _bam, bai -> [meta, bai] } // channel: [ val(meta), path(bai) ]
    aligned_reads_bam                   = val_skip_alignment ? channel.empty() : ch_aligned_bam.map { meta, bam, _bai -> [meta, bam] } // channel: [ val(meta), path(bam) ]
    aligned_reads_crai                  = !val_convert_unphased_aligned_reads_to_cram ? channel.empty() : SAMTOOLS_CONVERT.out.crai // channel: [ val(meta), path(crai) ]
    aligned_reads_cram                  = !val_convert_unphased_aligned_reads_to_cram ? channel.empty() : SAMTOOLS_CONVERT.out.cram // channel: [ val(meta), path(cram) ]
    assembly_summary                    = val_skip_genome_assembly ? channel.empty() : GENOME_ASSEMBLY.out.assembly_summary // channel: [ val(meta), path(assembly_summary) ]
    chromograph_plots                   = val_skip_chromograph ? channel.empty() : CHROMOGRAPH.out.chromograph_plots // channel: [ val(meta), path(png) ]
    cramino_phased_arrow                = val_skip_phasing ? channel.empty() : PHASING.out.haplotagging_arrow // channel: [ val(meta), path(arrow) ]
    cramino_phased_stats                = val_skip_phasing ? channel.empty() : PHASING.out.haplotagging_stats // channel: [ val(meta), path(txt) ]
    cramino_unphased_arrow              = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.cramino_arrow // channel: [ val(meta), path(arrow) ]
    cramino_unphased_stats              = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.cramino_stats // channel: [ val(meta), path(txt) ]
    fastqc_html                         = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.fastqc_html // channel: [ val(meta), path(html) ]
    fastqc_zip                          = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.fastqc_zip // channel: [ val(meta), path(zip) ]
    gens_baf_bed                        = val_skip_prepare_gens_input ? channel.empty() : PREPARE_GENS_INPUTS.out.baf_bed_tbi.map { meta, bed, _tbi -> [meta, bed] } // channel: [ val(meta), path(baf.bed.gz) ]
    gens_baf_tbi                        = val_skip_prepare_gens_input ? channel.empty() : PREPARE_GENS_INPUTS.out.baf_bed_tbi.map { meta, _bed, tbi -> [meta, tbi] } // channel: [ val(meta), path(baf.bed.gz.tbi) ]
    gens_cov_bed                        = val_skip_prepare_gens_input ? channel.empty() : PREPARE_GENS_INPUTS.out.cov_bed_tbi.map { meta, bed, _tbi -> [meta, bed] } // channel: [ val(meta), path(cov.bed.gz) ]
    gens_cov_tbi                        = val_skip_prepare_gens_input ? channel.empty() : PREPARE_GENS_INPUTS.out.cov_bed_tbi.map { meta, _bed, tbi -> [meta, tbi] } // channel: [ val(meta), path(cov.bed.gz.tbi) ]
    hificnv_copynum_bedgraph            = val_skip_sv_calling ? channel.empty() : CALL_SVS.out.hificnv_copynum // channel: [ val(meta), path(bedgraph) ]
    hificnv_depth_bw                    = val_skip_sv_calling ? channel.empty() : CALL_SVS.out.hificnv_depth // channel: [ val(meta), path(bw) ]
    hificnv_maf_bw                      = val_skip_sv_calling ? channel.empty() : CALL_SVS.out.hificnv_maf // channel: [ val(meta), path(bw) ]
    methylation_family_annotated        = val_skip_methylation_annotation ? channel.empty() : ANNOTATE_METHYLATION.out.methylation_annotation // channel: [ val(meta), path(methylated_regions_by_family) ]
    methylation_methbat_combined_bed    = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.pbcpg_combined_bed // channel: [ val(meta), path(bed.gz) ]
    methylation_methbat_combined_bigwig = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.pbcpg_combined_bigwig // channel: [ val(meta), path(combined.bw) ]
    methylation_methbat_combined_index  = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.pbcpg_combined_index // channel: [ val(meta), path(bed.gz.tbi) ]
    methylation_methbat_hap1_bed        = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.pbcpg_hap1_bed // channel: [ val(meta), path(bed.gz) ]
    methylation_methbat_hap1_bigwig     = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.pbcpg_hap1_bigwig // channel: [ val(meta), path(hap1.bw) ]
    methylation_methbat_hap1_index      = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.pbcpg_hap1_index // channel: [ val(meta), path(bed.gz.tbi) ]
    methylation_methbat_hap2_bed        = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.pbcpg_hap2_bed // channel: [ val(meta), path(bed.gz) ]
    methylation_methbat_hap2_bigwig     = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.pbcpg_hap2_bigwig // channel: [ val(meta), path(hap2.bw) ]
    methylation_methbat_hap2_index      = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.pbcpg_hap2_index // channel: [ val(meta), path(bed.gz.tbi) ]
    methylation_methbat_profiles        = val_skip_methbat ? channel.empty() : CALL_METHYLATION_METHBAT.out.region_profile // channel: [ val(meta), path(region_profile) ]
    methylation_modkit_bed              = val_skip_modkit ? channel.empty() : CALL_METHYLATION_MODKIT.out.bed // channel: [ val(meta), path(bed.gz) ]
    methylation_modkit_bigwig           = val_skip_modkit ? channel.empty() : CALL_METHYLATION_MODKIT.out.bigwig // channel: [ val(meta), path(bw) ]
    methylation_modkit_tbi              = val_skip_modkit ? channel.empty() : CALL_METHYLATION_MODKIT.out.tbi // channel: [ val(meta), path(bed.gz.tbi) ]
    mosdepth_global_dist                = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.mosdepth_global_dist // channel: [ val(meta), path(txt) ]
    mosdepth_per_base_d4                = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.mosdepth_per_base_d4 // channel: [ val(meta), path(d4) ]
    mosdepth_regions_bed                = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.mosdepth_regions_bed // channel: [ val(meta), path(bed.gz)     ]
    mosdepth_regions_csi                = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.mosdepth_regions_csi // channel: [ val(meta), path(bed.gz.csi) ]
    mosdepth_thresholds_bed             = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.mosdepth_thresholds_bed // channel: [ val(meta), path(bed.gz)     ]
    mosdepth_thresholds_csi             = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.mosdepth_thresholds_csi // channel: [ val(meta), path(bed.gz.csi) ]
    mosdepth_regions_dist               = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.mosdepth_regions_dist // channel: [ val(meta), path(txt) ]
    mosdepth_summary                    = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.mosdepth_summary // channel: [ val(meta), path(txt) ]
    multiqc_data                        = MULTIQC.out.data // channel: [ val(meta), path(*_data) ]
    multiqc_plots                       = MULTIQC.out.plots // channel: [ val(meta), path(*_plots) ]
    multiqc_report                      = MULTIQC.out.report // channel: [ val(meta), path(html) ]
    paralogs_family_annotated_json      = val_skip_annotate_paralogs ? channel.empty() : ANNOTATE_PARALOGS.out.json // channel: [ val(meta), path(json) ]
    paralogs_family_annotated_tsv       = val_skip_annotate_paralogs ? channel.empty() : ANNOTATE_PARALOGS.out.tsv // channel: [ val(meta), path(tsv) ]
    paralogs_family_tbi                 = val_skip_call_paralogs ? channel.empty() : CALL_PARALOGS.out.family_tbi // channel: [ val(meta), path(tbi) ]
    paralogs_family_vcf                 = val_skip_call_paralogs ? channel.empty() : CALL_PARALOGS.out.family_vcf // channel: [ val(meta), path(vcf) ]
    paralogs_sample_bai                 = val_skip_call_paralogs ? channel.empty() : CALL_PARALOGS.out.bai // channel: [ val(meta), path(bai) ]
    paralogs_sample_bam                 = val_skip_call_paralogs ? channel.empty() : CALL_PARALOGS.out.bam // channel: [ val(meta), path(bam) ]
    paralogs_sample_crai                = val_skip_call_paralogs ? channel.empty() : CALL_PARALOGS.out.crai // channel: [ val(meta), path(crai) ]
    paralogs_sample_cram                = val_skip_call_paralogs ? channel.empty() : CALL_PARALOGS.out.cram // channel: [ val(meta), path(cram) ]
    paralogs_sample_json                = val_skip_call_paralogs ? channel.empty() : CALL_PARALOGS.out.json // channel: [ val(meta), path(json) ]
    paralogs_sample_tbi                 = val_skip_call_paralogs ? channel.empty() : CALL_PARALOGS.out.sample_tbi // channel: [ val(meta), path(tbi) ]
    paralogs_sample_vcf                 = val_skip_call_paralogs ? channel.empty() : CALL_PARALOGS.out.sample_vcf // channel: [ val(meta), path(vcf) ]
    peddy_het_check_csv                 = val_skip_peddy ? channel.empty() : PEDDY.out.het_check_csv // channel: [ val(meta), path(csv) ]
    peddy_het_check_png                 = val_skip_peddy ? channel.empty() : PEDDY.out.het_check_png // channel: [ val(meta), path(png) ]
    peddy_html                          = val_skip_peddy ? channel.empty() : PEDDY.out.html // channel: [ val(meta), path(html) ]
    peddy_ped                           = val_skip_peddy ? channel.empty() : PEDDY.out.ped // channel: [ val(meta), path(ped) ]
    peddy_ped_check_csv                 = val_skip_peddy ? channel.empty() : PEDDY.out.ped_check_csv // channel: [ val(meta), path(csv) ]
    peddy_ped_check_png                 = val_skip_peddy ? channel.empty() : PEDDY.out.ped_check_png // channel: [ val(meta), path(png) ]
    peddy_ped_check_rel_difference_csv  = val_skip_peddy ? channel.empty() : PEDDY.out.ped_check_rel_difference_csv // channel: [ val(meta), path(csv) ]
    peddy_sex_check_csv                 = val_skip_peddy ? channel.empty() : PEDDY.out.sex_check_csv // channel: [ val(meta), path(csv) ]
    peddy_sex_check_png                 = val_skip_peddy ? channel.empty() : PEDDY.out.sex_check_png // channel: [ val(meta), path(png) ]
    peddy_vs_html                       = val_skip_peddy ? channel.empty() : PEDDY.out.vs_html // channel: [ val(meta), path(html) ]
    pedigree                            = val_skip_rank_variants ? channel.empty() : SOMALIER_PED_FAMILY.out.ped // channel: [ val(meta), path(ped) ]
    phasing_blocks_gtf                  = val_skip_phasing ? channel.empty() : PHASING.out.blocks // channel: [ val(meta), path("*.blocks.gtf.gz") ]
    phasing_blocks_tbi                  = val_skip_phasing ? channel.empty() : PHASING.out.blocks_index // channel: [ val(meta), path("*.blocks.gtf.gz.tbi") ]
    qc_bcftools_stats                   = val_skip_snv_calling ? channel.empty() : QC_SNVS.out.stats // channel: [ val(meta), path(txt) ]
    qc_deepvariant_vcfstatsreport       = val_skip_snv_calling ? channel.empty() : QC_SNVS.out.vcfstatsreport // channel: [ val(meta), path(html) ]
    qc_sambamba_depth_bed               = val_skip_qc ? channel.empty() : QC_ALIGNED_READS.out.sambamba_depth_bed // channel: [ val(meta), path(bed) ]
    qc_whatshap_stats                   = val_skip_phasing ? channel.empty() : PHASING.out.stats // channel: [ val(meta), path(*.stats.tsv) ]
    repeats_annotated_family_tbi        = val_skip_repeat_annotation ? channel.empty() : ANNOTATE_REPEAT_EXPANSIONS.out.tbi // channel: [ val(meta), path(tbi) ]
    repeats_annotated_family_vcf        = val_skip_repeat_annotation ? channel.empty() : ANNOTATE_REPEAT_EXPANSIONS.out.vcf // channel: [ val(meta), path(vcf) ]
    repeats_strdust_family_tbi          = val_skip_strdust ? channel.empty() : CALL_REPEAT_EXPANSIONS_STRDUST.out.family_tbi // channel: [ val(meta), path(tbi) ]
    repeats_strdust_family_vcf          = val_skip_strdust ? channel.empty() : CALL_REPEAT_EXPANSIONS_STRDUST.out.family_vcf // channel: [ val(meta), path(vcf) ]
    repeats_strdust_sample_tbi          = val_skip_strdust ? channel.empty() : CALL_REPEAT_EXPANSIONS_STRDUST.out.sample_tbi // channel: [ val(meta), path(tbi) ]
    repeats_strdust_sample_vcf          = val_skip_strdust ? channel.empty() : CALL_REPEAT_EXPANSIONS_STRDUST.out.sample_vcf // channel: [ val(meta), path(vcf) ]
    repeats_trgt_family_tbi             = val_skip_trgt ? channel.empty() : CALL_REPEAT_EXPANSIONS_TRGT.out.family_tbi // channel: [ val(meta), path(tbi) ]
    repeats_trgt_family_vcf             = val_skip_trgt ? channel.empty() : CALL_REPEAT_EXPANSIONS_TRGT.out.family_vcf // channel: [ val(meta), path(vcf) ]
    repeats_trgt_sample_bai             = val_skip_trgt ? channel.empty() : CALL_REPEAT_EXPANSIONS_TRGT.out.sample_bai // channel: [ val(meta), path(bai) ]
    repeats_trgt_sample_bam             = val_skip_trgt ? channel.empty() : CALL_REPEAT_EXPANSIONS_TRGT.out.sample_bam // channel: [ val(meta), path(bam) ]
    repeats_trgt_sample_crai            = val_skip_trgt ? channel.empty() : CALL_REPEAT_EXPANSIONS_TRGT.out.sample_crai // channel: [ val(meta), path(crai) ]
    repeats_trgt_sample_cram            = val_skip_trgt ? channel.empty() : CALL_REPEAT_EXPANSIONS_TRGT.out.sample_cram // channel: [ val(meta), path(cram) ]
    repeats_trgt_sample_tbi             = val_skip_trgt ? channel.empty() : CALL_REPEAT_EXPANSIONS_TRGT.out.sample_tbi // channel: [ val(meta), path(tbi) ]
    repeats_trgt_sample_vcf             = val_skip_trgt ? channel.empty() : CALL_REPEAT_EXPANSIONS_TRGT.out.sample_vcf // channel: [ val(meta), path(vcf) ]
    sawfish_copynum_bedgraph            = val_skip_sv_calling ? channel.empty() : CALL_SVS.out.sawfish_copynum_bedgraph // channel: [ val(meta), path(bedgraph) ]
    sawfish_depth_bw                    = val_skip_sv_calling ? channel.empty() : CALL_SVS.out.sawfish_depth_bw // channel: [ val(meta), path(bw) ]
    sawfish_gc_bias_corrected_depth_bw  = val_skip_sv_calling ? channel.empty() : CALL_SVS.out.sawfish_gc_bias_corrected_depth_bw // channel: [ val(meta), path(bw) ]
    sawfish_maf_bw                      = val_skip_sv_calling ? channel.empty() : CALL_SVS.out.sawfish_maf_bw // channel: [ val(meta), path(bw) ]
    somalier_relate_html                = val_skip_sex_check ? channel.empty() : BAM_INFER_SEX.out.somalier_html // channel: [ val(meta), path(html) ]
    somalier_relate_pairs               = val_skip_sex_check ? channel.empty() : BAM_INFER_SEX.out.somalier_pairs // channel: [ val(meta), path(pairs.tsv) ]
    somalier_relate_samples             = val_skip_sex_check ? channel.empty() : BAM_INFER_SEX.out.somalier_samples // channel: [ val(meta), path(samples.tsv) ]
    snvs_sample_tbi                     = val_skip_snv_calling ? channel.empty() : VCF_CONCAT_NORM_VARIANTS.out.index // channel: [ val(meta), path(tbi) ]
    snvs_sample_vcf                     = val_skip_snv_calling ? channel.empty() : VCF_CONCAT_NORM_VARIANTS.out.vcf // channel: [ val(meta), path(vcf) ]
    snvs_family_tbi                     = val_skip_snv_calling ? channel.empty() : CONCAT_SORT_RANKED_SNVS.out.index // channel: [ val(meta), path(tbi) ]
    snvs_family_vcf                     = val_skip_snv_calling ? channel.empty() : CONCAT_SORT_RANKED_SNVS.out.vcf // channel: [ val(meta), path(vcf) ]
    svs_per_family_and_caller_tbi       = val_skip_sv_calling ? channel.empty() : MERGE_SVS.out.family_caller_tbi // channel: [ val(meta), path(tbi) ]
    svs_per_family_and_caller_vcf       = val_skip_sv_calling ? channel.empty() : MERGE_SVS.out.family_caller_vcf // channel: [ val(meta), path(vcf) ]
    svs_per_family_tbi                  = val_skip_sv_calling ? channel.empty() : BCFTOOLS_VIEW_SV.out.tbi // channel: [ val(meta), path(tbi) ]
    svs_per_family_vcf                  = val_skip_sv_calling ? channel.empty() : BCFTOOLS_VIEW_SV.out.vcf // channel: [ val(meta), path(vcf.gz) ]
}

/**
 * Adds `child_with_two_parents_in_family` to the meta of `ch_input`, based on whether the family has a child with two parents according to the samplesheet.
 *
 * @param ch_input       Channel of [meta, file] where meta contains `family_id`
 * @param ch_samplesheet Channel of [meta, files] where meta contains `family_id` and `two_parents`
 * @return               Channel of [meta, file] with updated meta
 */
def addChildWithTwoParentsToMeta(ch_input, ch_samplesheet) {

    def ch_families = ch_samplesheet
        .map { meta, _files -> [meta.family_id, meta.two_parents] }
        .groupTuple()
        .map { family_id, child_with_two_parents -> [family_id, child_with_two_parents.any()] }

    // Need to use combine and filter here instead of join, since we only have one entry per family in ch_families but potentially multiple entries per family in ch_input
    ch_families
        .combine(ch_input)
        .filter { samplesheet_family_id, _child_with_two_parents, file_meta, _file ->
            samplesheet_family_id == file_meta.family_id
        }
        .map { _family_id, child_with_two_parents, meta, file ->
            [meta + [child_with_two_parents_in_family: child_with_two_parents], file]
        }
}

/**
 * Build input channel for ranking variants, by combining VCFs with PED files and ranking config, and adding variant type to the meta for downstream use.
 *
 * @param ch_vcf          Channel of [meta, vcf]
 * @param ch_ped          Channel of [meta, ped] (one per family)
 * @param variant_type    String (e.g. 'snv' or 'sv')
 * @param ch_score_config Channel of [meta, score_config]
 * @param ch_samplesheet  Channel used to derive family-level metadata
 * @return                Channel of [meta, vcf, ped, score_config]
 */
def buildRankVariantsInputChannel(ch_vcf, ch_ped, variant_type, ch_score_config, ch_samplesheet) {
    // This is used to determine compound ranking thresholds and penalties in genmod
    // Need to use combine and filter here instead of join, since we only have one entry per family in ch_ped but potentially multiple entries per family in vcf_with_meta.
    // The meta.id of the PED channel is the family_id
    addChildWithTwoParentsToMeta(ch_vcf, ch_samplesheet)
        .map { meta, vcf ->
            [meta + [variant_type: variant_type], vcf]
        }
        .combine(ch_ped)
        .filter { vcf_meta, _vcf, ped_meta, _ped ->
            vcf_meta.family_id == ped_meta.id
        }
        .combine(ch_score_config)
        .map { vcf_meta, vcf, _ped_meta, ped, _score_config_meta, score_config ->
            [vcf_meta, vcf, ped, score_config]
        }
}
