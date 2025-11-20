include { samplesheetToList } from 'plugin/nf-schema'
include {
    createReferenceChannelFromPath
    createReferenceChannelFromSamplesheet
} from '../subworkflows/local/utils_nfcore_nallo_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGN_ASSEMBLIES                                       } from '../subworkflows/local/align_assemblies'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV                    } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SVS                    } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_SNVS                                          } from '../subworkflows/local/annotate_snvs'
include { ANNOTATE_SVS                                           } from '../subworkflows/local/annotate_svs'
include { CONVERT_INPUT_FILES as CONVERT_INPUT_FASTQS            } from '../subworkflows/local/convert_input_files'
include { CONVERT_INPUT_FILES as CONVERT_INPUT_BAMS              } from '../subworkflows/local/convert_input_files'
include { BAM_INFER_SEX                                          } from '../subworkflows/local/bam_infer_sex'
include { CALL_PARALOGS                                          } from '../subworkflows/local/call_paralogs'
include { CALL_REPEAT_EXPANSIONS_STRDUST                         } from '../subworkflows/local/call_repeat_expansions_strdust'
include { CALL_REPEAT_EXPANSIONS_TRGT                            } from '../subworkflows/local/call_repeat_expansions_trgt'
include { CALL_SNVS                                              } from '../subworkflows/local/call_snvs'
include { CALL_SVS                                               } from '../subworkflows/local/call_svs'
include { GENOME_ASSEMBLY                                        } from '../subworkflows/local/genome_assembly'
include { GVCF_GLNEXUS_NORM_VARIANTS                             } from '../subworkflows/local/gvcf_glnexus_norm_variants'
include { METHYLATION                                            } from '../subworkflows/local/methylation'
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
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// local
include { CREATE_PEDIGREE_FILE as SAMPLESHEET_PED                } from '../modules/local/create_pedigree_file/main'
include { CREATE_PEDIGREE_FILE as SOMALIER_PED                   } from '../modules/local/create_pedigree_file/main'
include { CREATE_PEDIGREE_FILE as SOMALIER_PED_FAMILY            } from '../modules/local/create_pedigree_file/main'

// nf-core
include { BCFTOOLS_CONCAT                                        } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT                                          } from '../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_VIEW                                          } from '../modules/nf-core/bcftools/view/main'
include { MINIMAP2_ALIGN                                         } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                                         } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_CONVERT                                       } from '../modules/nf-core/samtools/convert/main'
include { MULTIQC                                                } from '../modules/nf-core/multiqc/main'
include { PEDDY                                                  } from '../modules/nf-core/peddy/main'
include { SPLITUBAM                                              } from '../modules/nf-core/splitubam/main'
include { STRANGER                                               } from '../modules/nf-core/stranger/main'
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
    ch_input

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Channels from (optional) input files
    // If provided: [[id: 'reference'], [/path/to/reference_full_name.file]]
    ch_cadd_header               = createReferenceChannelFromPath("$projectDir/assets/cadd_to_vcf_header_-1.0-.txt")
    ch_cadd_resources            = createReferenceChannelFromPath(params.cadd_resources)
    ch_cadd_prescored_indels     = createReferenceChannelFromPath(params.cadd_prescored_indels)
    ch_fasta                     = createReferenceChannelFromPath(params.fasta)
    ch_tandem_repeats            = createReferenceChannelFromPath(params.tandem_repeats, Channel.value([[],[]]))
    ch_par                       = createReferenceChannelFromPath(params.par_regions)
    ch_str_bed                   = createReferenceChannelFromPath(params.str_bed)
    ch_snv_call_regions          = createReferenceChannelFromPath(params.snv_call_regions, Channel.value([[],[]]))
    ch_sv_call_regions           = createReferenceChannelFromPath(params.sv_call_regions)
    ch_methylation_call_regions  = createReferenceChannelFromPath(params.methylation_call_regions, Channel.value([[],[]]))
    ch_stranger_repeat_catalog   = createReferenceChannelFromPath(params.stranger_repeat_catalog)
    ch_variant_consequences_snvs = createReferenceChannelFromPath(params.variant_consequences_snvs)
    ch_variant_consequences_svs  = createReferenceChannelFromPath(params.variant_consequences_svs)
    ch_vep_cache_unprocessed     = createReferenceChannelFromPath(params.vep_cache, Channel.value([[],[]]))
    ch_expected_xy_bed           = createReferenceChannelFromPath(params.cnv_expected_xy_cn)
    ch_expected_xx_bed           = createReferenceChannelFromPath(params.cnv_expected_xx_cn)
    ch_exclude_bed               = createReferenceChannelFromPath(params.cnv_excluded_regions)
    ch_genmod_reduced_penetrance = createReferenceChannelFromPath(params.genmod_reduced_penetrance)
    ch_genmod_score_config_snvs  = createReferenceChannelFromPath(params.genmod_score_config_snvs)
    ch_genmod_score_config_svs   = createReferenceChannelFromPath(params.genmod_score_config_svs)
    ch_peddy_sites               = createReferenceChannelFromPath(params.peddy_sites, Channel.value([[],[]]))
    ch_qc_regions                = createReferenceChannelFromPath(params.qc_regions, Channel.value([[],[]]))
    ch_somalier_sites            = createReferenceChannelFromPath(params.somalier_sites)


    // Channels from (optional) input samplesheets validated by schema
    ch_databases                 = createReferenceChannelFromSamplesheet(params.echtvar_snv_databases, 'assets/schema_snp_db.json', Channel.value([[],[]]))
    ch_svdb_sv_databases         = createReferenceChannelFromSamplesheet(params.svdb_sv_databases, 'assets/svdb_query_vcf_schema.json', Channel.value([]))
    ch_vep_plugin_files          = createReferenceChannelFromSamplesheet(params.vep_plugin_files, 'assets/schema_vep_plugin_files.json', Channel.value([]))
    ch_hgnc_ids                  = createReferenceChannelFromSamplesheet(params.filter_variants_hgnc_ids, 'assets/schema_hgnc_ids.json', Channel.value([]))
        .map { it[0].toString() } // only one element per row
        .collectFile(name: 'hgnc_ids.txt', newLine: true, sort: true)
        .map { file -> [ [ id: 'hgnc_ids' ], file ] }
        .collect()

    def cram_output = params.alignment_output_format == 'cram'

    //
    // Prepare references
    //
    if(!params.skip_alignment || !params.skip_genome_assembly) {

        // The genome assembly alignment needs a fai for cram output,
        // but we shouldn't need to prepare the vep cache here.
        // Perhaps PREPARE_REFERENCES could be modified to handle this case?
        PREPARE_REFERENCES (
            ch_fasta,
            ch_vep_cache_unprocessed,
            params.fasta.endsWith('.gz'),                           // should we unzip fasta
            params.vep_cache && params.vep_cache.endsWith("tar.gz") // should we untar vep cache
        )
        ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

        // Gather indices
        ch_fasta = PREPARE_REFERENCES.out.fasta
        ch_fai   = PREPARE_REFERENCES.out.fai
    }

    // Convert FASTQ to BAM only if alignment or should be done.
    // Since we assume that the majority of pipeline runs will use BAM files as input,
    // we start all files as BAMs for simplicity except for the assembly, which requires FASTQs.
    if(!params.skip_alignment) {

        CONVERT_INPUT_FASTQS (
            ch_input,
            false,
            true
        )
        ch_versions = ch_versions.mix(CONVERT_INPUT_FASTQS.out.versions)
    }

    // To speed up the alignement, we can split the BAM files into smaller chunks.
    // We can also use the split BAM files for FASTQ conversion to the assembly workflow,
    // instead of the original BAM files which should allow the assembly to start sooner.
    //
    // We could change the name of alignment processes to something more generic, like `--split_input_files`?
    // If we start running more trios we also need to consider that the parents at the moment needs to be merged
    // before YAK. So, we could consider adding some logic to handle that case,
    // to avoid unneccessary splitting and merging just for a minor speedup in the conversion.
    if(!params.skip_alignment && params.alignment_processes > 1) {

        SPLITUBAM (
            CONVERT_INPUT_FASTQS.out.bam // contains all BAM files, including those not converted.
        )
        ch_versions = ch_versions.mix(SPLITUBAM.out.versions)

    }

    //
    // Hifiasm assembly and alignment to reference genome
    //
    if(!params.skip_genome_assembly) {

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

        CONVERT_INPUT_BAMS (
            !params.skip_alignment && params.alignment_processes > 1 ? SPLITUBAM.out.bam.transpose() : ch_input,
            true,
            false
        )
        ch_versions = ch_versions.mix(CONVERT_INPUT_BAMS.out.versions)

        // contains all FASTQ files, including those not converted
        CONVERT_INPUT_BAMS.out.fastq
            .groupTuple()
            .set { ch_genome_assembly_input }

        // Hifiasm assembly
        GENOME_ASSEMBLY(
            ch_genome_assembly_input,
            params.hifiasm_mode == "trio-binning" // Should we use trio binning mode?
        )
        ch_versions = ch_versions.mix(GENOME_ASSEMBLY.out.versions)

        ALIGN_ASSEMBLIES (
            GENOME_ASSEMBLY.out.assembled_haplotypes,
            ch_fasta,
            ch_fai,
            cram_output
        )
        ch_versions = ch_versions.mix(ALIGN_ASSEMBLIES.out.versions)
    }

    //
    // Map reads to reference
    //
    if (!params.skip_alignment) {

        (params.alignment_processes > 1 ? SPLITUBAM.out.bam.transpose() : CONVERT_INPUT_FASTQS.out.bam)
            .set { reads_for_alignment }

        // If we are splitting input files, this needs to wait for all files to be split. But while we do that, the reads can still start to be aligned.
        // Later, this allows us to combine this meta with the aligned BAMs for merging, without having to wait for all alignment processes to finish.
        reads_for_alignment
            .groupTuple()
            .map { meta, files -> [ meta + [ n_files: files.size() ] ] }
            .set { reads_grouping_key }

        //
        // Align reads (could be a split-align-merge subworkflow)
        //
        MINIMAP2_ALIGN (
            reads_for_alignment,
            PREPARE_REFERENCES.out.mmi,
            true,
            'bai',
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        // Split channel into cases where we have multiple files or single files
        MINIMAP2_ALIGN.out.bam
            // If there are multiple files per sample, each file has the same meta so failOnDuplicate fails here.
            // The end result is fine, but it might be worth to e.g. give each file a non-identical meta,
            // then join, strip identifier, join again, to be able to run the pipeline in strict mode.
            .join(MINIMAP2_ALIGN.out.index, failOnMismatch:true)
            .combine(reads_grouping_key)
            .filter { aligned_meta, _bam, _bai, grouping_key_meta ->
                aligned_meta.id == grouping_key_meta.id
            }
            .map { aligned_meta, bam, bai, grouping_key_meta ->
                [ aligned_meta + [ n_files: grouping_key_meta.n_files ], bam, bai ]
            }
            .map { meta, bam, bai ->
                [ groupKey(meta, meta.n_files), bam, bai ]
            }
            .groupTuple()
            .set { bam_to_merge }

        // Merge files - even if we only have one file per sample.
        // This is because sometimes we need to output unphased BAM files. We can no longer output single files
        // from the alignment process, because they would need to be renamed, based on n_files, which is no longer
        // available to the alignment process. Because that would mean having to wait for all samples to be alinged
        // before moving on to subsequent steps.
        SAMTOOLS_MERGE (
            bam_to_merge.map { meta, bam, _bai -> [ meta, bam ] },
            [[],[]],
            [[],[]],
            [[],[]],
        )
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

        // Combine merged with unmerged bam files
        SAMTOOLS_MERGE.out.bam
            .join(SAMTOOLS_MERGE.out.bai, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, bam, bai -> [ meta - meta.subMap('n_files'), bam, bai ] }
            .set { ch_aligned_bam }

        // Publish alignments as CRAM if requested
        if (cram_output && params.skip_phasing) {
            SAMTOOLS_CONVERT (
                ch_aligned_bam,
                ch_fasta,
                ch_fai
            )
            ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)
        }

        //
        // Create PED from samplesheet
        //
        ch_input
            .map { meta, _files -> [ [ id: meta.project ], meta ] }
            .groupTuple()
            .set { ch_samplesheet_ped_in }

        SAMPLESHEET_PED ( ch_samplesheet_ped_in )
        ch_versions = ch_versions.mix(SAMPLESHEET_PED.out.versions)

        SAMPLESHEET_PED.out.ped
            .collect()
            .set { ch_samplesheet_pedfile }

        //
        // Check sex and relatedness, and update with inferred sex if the sex for a sample is unknown
        //
        BAM_INFER_SEX (
            ch_aligned_bam,
            ch_fasta,
            ch_fai,
            ch_somalier_sites,
            ch_samplesheet_pedfile
        )
        ch_versions = ch_versions.mix(BAM_INFER_SEX.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_INFER_SEX.out.somalier_samples.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BAM_INFER_SEX.out.somalier_pairs.map{it[1]}.collect().ifEmpty([]))

        // Set files with updated meta for subsequent processes
        ch_bam     = BAM_INFER_SEX.out.bam
        ch_bam_bai = BAM_INFER_SEX.out.bam_bai

    }

    //
    // Run read QC with FastQC, mosdepth and cramino
    //
    if (!params.skip_qc) {

        QC_ALIGNED_READS (
            ch_bam_bai,
            ch_fasta,
            ch_qc_regions,
        )
        ch_versions = ch_versions.mix(QC_ALIGNED_READS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix( QC_ALIGNED_READS.out.fastqc_zip.collect { it[1] }.ifEmpty([]) )
        ch_multiqc_files = ch_multiqc_files.mix( QC_ALIGNED_READS.out.mosdepth_summary.collect { it[1] } )
        ch_multiqc_files = ch_multiqc_files.mix( QC_ALIGNED_READS.out.mosdepth_global_dist.collect { it[1] } )
        ch_multiqc_files = ch_multiqc_files.mix( QC_ALIGNED_READS.out.mosdepth_region_dist.collect { it[1] }.ifEmpty([]) )

    }

    //
    // Call paralogous genes with paraphase
    //
    if(!params.skip_call_paralogs) {
        CALL_PARALOGS (
            ch_bam_bai,
            ch_fasta,
            ch_fai,
            cram_output
        )
        ch_versions = ch_versions.mix(CALL_PARALOGS.out.versions)
    }

    //
    // Call SNVs
    //
    if(!params.skip_snv_calling) {

        // Make BED intervals, can be used for parallel SNV calling
        SCATTER_GENOME (
            ch_fai,
            ch_snv_call_regions,      // BED file to scatter
            !params.snv_call_regions, // Make bed from fai
            !params.skip_snv_calling,
            params.snv_calling_processes
        )
        ch_versions = ch_versions.mix(SCATTER_GENOME.out.versions)

        // Combine the BED intervals with BAM/BAI files to create a region-bam-bai for each sample.
        // This uses the whole BAM files for each region instead of splitting them.
        ch_bam_bai
            .combine(SCATTER_GENOME.out.bed_intervals)
            .map { meta, bam, bai, bed, intervals ->
                [ meta + [ num_intervals: intervals, region: bed ], bam, bai, bed ]
            }
            .set { call_snvs_input }

        CALL_SNVS(
            call_snvs_input,
            ch_fasta,
            ch_fai,
            ch_par,
            "deepvariant",
        )
        ch_versions = ch_versions.mix(CALL_SNVS.out.versions)

        CALL_SNVS.out.gvcf
            .map { meta, gvcf ->
                [[id: meta.region.name, family_id: meta.family_id], gvcf]
            }
            .groupTuple()
            .set { variants_to_merge_per_family }

        // Create a merged and normalized VCF, containing one region with all samples, to be used in annotation and ranking.
        GVCF_GLNEXUS_NORM_VARIANTS(
            variants_to_merge_per_family,
            SCATTER_GENOME.out.bed, // This contains all regions, but we could probably pass the region BED that actually matches the variants instead...
            ch_fasta,
            "deepvariant",
        )
        ch_versions = ch_versions.mix(GVCF_GLNEXUS_NORM_VARIANTS.out.versions)

        if (params.prepare_gens_input) {

            def missingGensParams = [
                params.gens_baf_positions ? null : '--gens_baf_positions',
                params.gens_gatk_header_template ? null : '--gens_gatk_header_template',
                params.gens_panel_of_normals ? null : '--gens_panel_of_normals',
            ].findAll { it }

            if (missingGensParams) {
                error "The following parameters are required when --prepare_gens_input is enabled: ${missingGensParams.join(', ')}"
            }

            CALL_SNVS.out.gvcf
                .join(CALL_SNVS.out.gvcf_index)
                .map { meta, gvcf, gvcf_index -> 
                    def sample_meta = meta - meta.subMap('region')
                    [ sample_meta, gvcf, gvcf_index ]
                }
                .groupTuple()
                .set { ch_gvcfs_to_concat_per_sample }

            // ch_gvcfs_to_concat_per_sample.view()

            def ch_gens_baf_positions = Channel.fromPath(params.gens_baf_positions, checkIfExists: true)
            def ch_gens_gatk_header_template = Channel.fromPath(params.gens_gatk_header_template, checkIfExists: true)
            def ch_gens_panel_of_normals = Channel.fromPath(params.gens_panel_of_normals, checkIfExists: true)

            PREPARE_GENS_INPUTS(
               ch_bam_bai,
               ch_gvcfs_to_concat_per_sample,
               ch_gens_baf_positions,
               ch_gens_gatk_header_template,
               ch_gens_panel_of_normals
            )
            ch_versions = ch_versions.mix(PREPARE_GENS_INPUTS.out.versions)
        }

        CALL_SNVS.out.vcf
            .map { meta, vcf ->
                def new_meta = meta - meta.subMap('region')
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
            "deepvariant",
        )
        ch_versions = ch_versions.mix(VCF_CONCAT_NORM_VARIANTS.out.versions)

        // These contains RefCalls
        sample_snv_vcf   = VCF_CONCAT_NORM_VARIANTS.out.vcf
        sample_snv_index = VCF_CONCAT_NORM_VARIANTS.out.index

        family_snv_vcf   = GVCF_GLNEXUS_NORM_VARIANTS.out.vcf
        family_snv_index = GVCF_GLNEXUS_NORM_VARIANTS.out.index

        // SNV QC
        QC_SNVS (
            VCF_CONCAT_NORM_VARIANTS.out.bcftools_concat_vcf, // Can we use the normalized VCF here, for DV vcfstatsreport?
            sample_snv_vcf,
            sample_snv_index,
        )

        ch_versions = ch_versions.mix(QC_SNVS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(QC_SNVS.out.stats.collect{it[1]}.ifEmpty([]))

        family_snv_vcf
            .join(family_snv_index, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_vcf_tbi_per_region }
    }

    //
    // Annotate SNVs
    //
    if(!params.skip_snv_annotation) {

        // Annotates family VCFs per variant call region
        ANNOTATE_SNVS(
            family_snv_vcf,
            ch_databases.map { _meta, databases -> databases }.collect(),
            ch_fasta,
            ch_fai,
            PREPARE_REFERENCES.out.vep_resources.map { _meta, cache -> cache },
            params.vep_cache_version,
            ch_vep_plugin_files.collect(),
            params.cadd_resources && params.cadd_prescored_indels,
            params.echtvar_snv_databases,
            ch_cadd_header,
            ch_cadd_resources,
            ch_cadd_prescored_indels,
            params.pre_vep_snv_filter_expression != ''
        )
        ch_versions = ch_versions.mix(ANNOTATE_SNVS.out.versions)

        ANNOTATE_SNVS.out.vcf
            .multiMap { meta, vcf ->
                clinical: [ meta + [ set: "clinical" ], vcf ]
                research: [ meta + [ set: "research" ], vcf ]
            }
            .set { ch_clin_research_snvs_vcf }

        ch_clin_research_snvs_vcf.research
            .set { ch_ann_csq_pli_snv_in }

        if(params.filter_variants_hgnc_ids || params.filter_snvs_expression != '') {

            FILTER_VARIANTS_SNVS (
                ch_clin_research_snvs_vcf.clinical,
                ch_hgnc_ids,
                params.filter_snvs_expression,
                params.filter_variants_hgnc_ids,
            )
            ch_versions = ch_versions.mix(FILTER_VARIANTS_SNVS.out.versions)

            ch_ann_csq_pli_snv_in = ch_ann_csq_pli_snv_in.mix(FILTER_VARIANTS_SNVS.out.vcf)
        }

        ANN_CSQ_PLI_SNV (
            ch_ann_csq_pli_snv_in,
            ch_variant_consequences_snvs
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SNV.out.versions)

        ANN_CSQ_PLI_SNV.out.vcf
            .join( ANN_CSQ_PLI_SNV.out.tbi, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_vcf_tbi_per_region }

    }

    //
    // Ranks family VCFs per variant call region
    // Can only run if samplesheet has affected samples
    //
    if(!params.skip_rank_variants) {

        // Create PED with updated sex - per family
        SOMALIER_PED_FAMILY (
            ch_bam
                .map { meta, _files -> [ [ id: meta.family_id ], meta ] }
                .groupTuple()
        )
        ch_versions = ch_versions.mix(SOMALIER_PED_FAMILY.out.versions)

        // Give PED file SNV meta so they can be joined later in the subworkflow.
        // Since we don't always have matching number of ped files and call regions
        // we need to combine and filter instead of join
        ANN_CSQ_PLI_SNV.out.vcf
            .map { meta, _vcf -> [ [ id:meta.family_id ], meta ] }
            .combine ( SOMALIER_PED_FAMILY.out.ped )
            .filter { family_id_snv, _meta, family_id_ped, _ped -> family_id_snv == family_id_ped }
            .map { _family_id_snv, meta, _family_id_ped, ped -> [ meta, ped ] }
            .set { ch_snv_ranking_ped_file }

        // Only run if we have affected individuals
        RANK_VARIANTS_SNV (
            ANN_CSQ_PLI_SNV.out.vcf,
            ch_snv_ranking_ped_file,
            ch_genmod_reduced_penetrance,
            ch_genmod_score_config_snvs
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)

        RANK_VARIANTS_SNV.out.vcf
            .join( RANK_VARIANTS_SNV.out.tbi, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_vcf_tbi_per_region }
    }

    //
    // Concatenate and sort SNVs (could be a subworkflow)
    //
    if(!params.skip_snv_calling) {

        ch_vcf_tbi_per_region
            .map { meta, vcf, tbi -> [ [ id: meta.family_id, set: meta.set ], vcf, tbi ] }
            .groupTuple(size: params.snv_calling_processes)
            .set { ch_bcftools_concat_in }

        // Concat into family VCFs per family with all regions
        BCFTOOLS_CONCAT (
                ch_bcftools_concat_in
            )
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

        // Sort and publish
        BCFTOOLS_SORT ( BCFTOOLS_CONCAT.out.vcf )
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)
    }

    //
    // Run Peddy
    //
    if (!params.skip_snv_calling && !params.skip_peddy) {

        BCFTOOLS_SORT.out.vcf
            .join( BCFTOOLS_SORT.out.tbi )
            .filter { meta, _vcf, _tbi -> meta.set == "research" }
            .set { ch_peddy_in }

        PEDDY (
            ch_peddy_in,
            ch_samplesheet_pedfile,
            ch_peddy_sites
        )
        ch_versions = ch_versions.mix(PEDDY.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.het_check_csv.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.sex_check_csv.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped_check_csv.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped_check_rel_difference_csv.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.het_check_png.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.sex_check_png.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped_check_png.map{it[1]}.collect().ifEmpty([]))

    }

    //
    // Call SVs
    //
    if(!params.skip_sv_calling) {

        CALL_SVS (
            ch_bam_bai,
            ch_tandem_repeats,
            sample_snv_vcf,
            ch_fasta,
            ch_expected_xy_bed,
            ch_expected_xx_bed,
            ch_exclude_bed,
            params.sv_callers_to_run.split(',').collect { it.toLowerCase().trim() },
            params.sv_callers_to_merge.split(',').collect { it.toLowerCase().trim() },
            params.sv_callers_merge_priority.split(',').collect { it.toLowerCase().trim() },
            ch_sv_call_regions,
            params.sv_call_regions,
            params.force_sawfish_joint_call_single_samples,
        )

        ch_versions = ch_versions.mix(CALL_SVS.out.versions)

    }

    //
    // Annotate SVs
    //
    if (!params.skip_sv_annotation) {

        ANNOTATE_SVS (
            CALL_SVS.out.family_vcf,
            ch_fasta,
            ch_svdb_sv_databases,
            PREPARE_REFERENCES.out.vep_resources.map { _meta, cache -> cache },
            params.vep_cache_version,
            ch_vep_plugin_files.collect()
        )
        ch_versions = ch_versions.mix(ANNOTATE_SVS.out.versions)

        ANNOTATE_SVS.out.vcf
            .multiMap { meta, vcf ->
                clinical: [ meta + [ set: "clinical" ], vcf ]
                research: [ meta + [ set: "research" ], vcf ]
            }
            .set { ch_clin_research_svs_vcf }

        ch_clin_research_svs_vcf.research
            .set { ch_ann_csq_svs_in }

        //
        // Filter SVs
        //
        if(params.filter_variants_hgnc_ids || params.filter_svs_expression != '') {

            FILTER_VARIANTS_SVS (
                ch_clin_research_svs_vcf.clinical,
                ch_hgnc_ids,
                params.filter_svs_expression,
                params.filter_variants_hgnc_ids,
            )
            ch_versions = ch_versions.mix(FILTER_VARIANTS_SVS.out.versions)

            ch_ann_csq_svs_in = ch_ann_csq_svs_in.mix(FILTER_VARIANTS_SVS.out.vcf)

        }

        ANN_CSQ_PLI_SVS (
            ch_ann_csq_svs_in,
            ch_variant_consequences_svs
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SVS.out.versions)
    }

    //
    // Rank SVs
    //
    if (!params.skip_rank_variants) {

        // Give PED file SVs meta so they can be joined later in the subworkflow.
        ANN_CSQ_PLI_SVS.out.vcf
            .combine ( SOMALIER_PED_FAMILY.out.ped )
            .filter { vcf_meta, _vcf, ped_meta, _ped -> vcf_meta.id == ped_meta.id }
            .map { vcf_meta, _vcf, _ped_meta, ped -> [ vcf_meta, ped ] }
            .set { ch_sv_ranking_ped_file }

        RANK_VARIANTS_SVS (
            ANN_CSQ_PLI_SVS.out.vcf,
            ch_sv_ranking_ped_file,
            ch_genmod_reduced_penetrance,
            ch_genmod_score_config_svs
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SVS.out.versions)
    }

    //
    // Collect and publish SVs
    //
    if(!params.skip_sv_calling) {

        ch_collect_svs = params.skip_sv_annotation ? CALL_SVS.out.family_vcf :
            params.skip_rank_variants ? ANN_CSQ_PLI_SVS.out.vcf :
            RANK_VARIANTS_SVS.out.vcf

        BCFTOOLS_VIEW (
            ch_collect_svs.map { meta, vcf -> [ meta, vcf, [] ] },
            [],
            [],
            []
        )
    }

    //
    // Phase SNVs and INDELs
    //
    if(!params.skip_phasing) {

        PHASING (
            sample_snv_vcf,
            sample_snv_index,
            ch_bam_bai,
            ch_fasta,
            ch_fai,
            cram_output
        )
        ch_versions = ch_versions.mix(PHASING.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(PHASING.out.stats.collect{it[1]}.ifEmpty([]))

    }

    //
    // Create methylation pileups with modkit
    //
    if(!params.skip_methylation_pileups) {
        METHYLATION (
            !params.skip_phasing ? PHASING.out.haplotagged_bam_bai : ch_bam_bai,
            ch_fasta,
            ch_fai,
            ch_methylation_call_regions,
            params.bigwig_modcodes
        )
        ch_versions = ch_versions.mix(METHYLATION.out.versions)
    }

    //
    // Call repeat expansions with TRGT
    //
    if(!params.skip_repeat_calling) {
        if (params.str_caller == "trgt") {
            CALL_REPEAT_EXPANSIONS_TRGT (
                PHASING.out.haplotagged_bam_bai,
                ch_fasta,
                ch_fai,
                ch_str_bed,
                cram_output
            )
            ch_versions = ch_versions.mix(CALL_REPEAT_EXPANSIONS_TRGT.out.versions)
            ch_repeat_expansions = CALL_REPEAT_EXPANSIONS_TRGT.out.family_vcf
        } else if (params.str_caller == "strdust"){
            CALL_REPEAT_EXPANSIONS_STRDUST (
                PHASING.out.haplotagged_bam_bai,
                ch_fasta,
                ch_fai,
                ch_str_bed
            )
            ch_versions = ch_versions.mix(CALL_REPEAT_EXPANSIONS_STRDUST.out.versions)
        }
    }

    //
    // Annotate repeat expansions with Stranger
    //
    if(!params.skip_repeat_annotation) {
        STRANGER (
            ch_repeat_expansions,
            ch_stranger_repeat_catalog
        )
        ch_versions = ch_versions.mix(STRANGER.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nallo_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.of(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )
    ch_methods_description_citation       = citationBibliographyText(
        ch_versions,
        file("$projectDir/assets/software_references.yml"),
        'citation'
    )
    ch_methods_description_bibliography   = citationBibliographyText(
        ch_versions,
        file("$projectDir/assets/software_references.yml"),
        'bibliography'
    )
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description
            .concat(ch_methods_description_citation)
            .concat(ch_methods_description_bibliography)
            .flatten()
            .collectFile(
                name: 'methods_description_mqc.yaml',
                sort: false // preserve order for correct yaml structure
            )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
