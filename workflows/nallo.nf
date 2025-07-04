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

include { ALIGN_ASSEMBLIES                        } from '../subworkflows/local/align_assemblies'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV     } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SVS     } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_SVS                            } from '../subworkflows/local/annotate_svs'
include { ASSEMBLY                                } from '../subworkflows/local/genome_assembly'
include { CONVERT_INPUT_FILES                     } from '../subworkflows/local/convert_input_files'
include { BAM_INFER_SEX                           } from '../subworkflows/local/bam_infer_sex'
include { CALL_CNVS                               } from '../subworkflows/local/call_cnvs'
include { CALL_PARALOGS                           } from '../subworkflows/local/call_paralogs'
include { CALL_REPEAT_EXPANSIONS_STRDUST          } from '../subworkflows/local/call_repeat_expansions_strdust'
include { CALL_REPEAT_EXPANSIONS_TRGT             } from '../subworkflows/local/call_repeat_expansions_trgt'
include { CALL_SVS                                } from '../subworkflows/local/call_svs'
include { FILTER_VARIANTS as FILTER_VARIANTS_SNVS } from '../subworkflows/local/filter_variants'
include { FILTER_VARIANTS as FILTER_VARIANTS_SVS  } from '../subworkflows/local/filter_variants'
include { METHYLATION                             } from '../subworkflows/local/methylation'
include { PHASING                                 } from '../subworkflows/local/phasing'
include { PREPARE_GENOME                          } from '../subworkflows/local/prepare_genome'
include { QC_ALIGNED_READS                        } from '../subworkflows/local/qc_aligned_reads'
include { RANK_VARIANTS as RANK_VARIANTS_SNV      } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SVS      } from '../subworkflows/local/rank_variants'
include { SCATTER_GENOME                          } from '../subworkflows/local/scatter_genome'
include { SHORT_VARIANT_CALLING                   } from '../subworkflows/local/short_variant_calling'
include { SNV_ANNOTATION                          } from '../subworkflows/local/snv_annotation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// local
include { CREATE_PEDIGREE_FILE as SAMPLESHEET_PED           } from '../modules/local/create_pedigree_file/main'
include { CREATE_PEDIGREE_FILE as SOMALIER_PED              } from '../modules/local/create_pedigree_file/main'
include { CREATE_PEDIGREE_FILE as SOMALIER_PED_FAMILY       } from '../modules/local/create_pedigree_file/main'

// nf-core
include { BCFTOOLS_CONCAT                                   } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT                                     } from '../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_STATS                                    } from '../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_VIEW                                     } from '../modules/nf-core/bcftools/view/main'
include { MINIMAP2_ALIGN                                    } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                                    } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_CONVERT                                  } from '../modules/nf-core/samtools/convert/main'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { PEDDY                                             } from '../modules/nf-core/peddy/main'
include { SPLITUBAM                                         } from '../modules/nf-core/splitubam/main'
include { STRANGER                                          } from '../modules/nf-core/stranger/main'
include { SVDB_MERGE as SVDB_MERGE_SVS_CNVS                 } from '../modules/nf-core/svdb/merge/main'
include { paramsSummaryMap                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../subworkflows/local/utils_nfcore_nallo_pipeline'
include { citationBibliographyText                          } from '../subworkflows/local/utils_nfcore_nallo_pipeline'

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
    ch_input_bed                 = createReferenceChannelFromPath(params.target_regions, Channel.value([[],[]]))
    ch_par                       = createReferenceChannelFromPath(params.par_regions)
    ch_str_bed                   = createReferenceChannelFromPath(params.str_bed)
    ch_stranger_repeat_catalog   = createReferenceChannelFromPath(params.stranger_repeat_catalog)
    ch_variant_consequences_snvs = createReferenceChannelFromPath(params.variant_consequences_snvs)
    ch_variant_consequences_svs  = createReferenceChannelFromPath(params.variant_consequences_svs)
    ch_vep_cache_unprocessed     = createReferenceChannelFromPath(params.vep_cache, Channel.value([[],[]]))
    ch_expected_xy_bed           = createReferenceChannelFromPath(params.hificnv_expected_xy_cn)
    ch_expected_xx_bed           = createReferenceChannelFromPath(params.hificnv_expected_xx_cn)
    ch_exclude_bed               = createReferenceChannelFromPath(params.hificnv_excluded_regions)
    ch_genmod_reduced_penetrance = createReferenceChannelFromPath(params.genmod_reduced_penetrance)
    ch_genmod_score_config_snvs  = createReferenceChannelFromPath(params.genmod_score_config_snvs)
    ch_genmod_score_config_svs   = createReferenceChannelFromPath(params.genmod_score_config_svs)
    ch_peddy_sites               = createReferenceChannelFromPath(params.peddy_sites, Channel.value([[],[]]))
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
    // Convert FASTQ to BAM (and vice versa if assembly workflow is active)
    //
    CONVERT_INPUT_FILES (
        ch_input,
        !params.skip_genome_assembly, // should bam -> fastq conversion be done
        !params.skip_alignment        // should fastq -> bam conversion be done
    )
    ch_versions = ch_versions.mix(CONVERT_INPUT_FILES.out.versions)

    //
    // Map reads to reference
    //
    if (!params.skip_alignment) {

        // Prepare references
        PREPARE_GENOME (
            ch_fasta,
            ch_vep_cache_unprocessed,
            params.fasta.endsWith('.gz'),                           // should we unzip fasta
            params.vep_cache && params.vep_cache.endsWith("tar.gz") // should we untar vep cache
        )
        ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

        // Gather indices
        ch_fasta = PREPARE_GENOME.out.fasta
        ch_fai   = PREPARE_GENOME.out.fai

        // Split input files for alignment
        if (params.alignment_processes > 1) {

            SPLITUBAM ( CONVERT_INPUT_FILES.out.bam )
            ch_versions = ch_versions.mix(SPLITUBAM.out.versions)
        }

        //
        // Align reads (could be a split-align-merge subworkflow)
        //
        MINIMAP2_ALIGN (
            params.alignment_processes > 1 ? SPLITUBAM.out.bam.transpose() : CONVERT_INPUT_FILES.out.bam,
            PREPARE_GENOME.out.mmi,
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
            .map {
                meta, bam, bai ->
                    [ groupKey(meta, meta.n_files), bam, bai ]
            }
            .groupTuple()
            .branch { meta, bam, bai ->
                single:   meta.n_files <= 1
                    return [ meta, bam[0], bai[0] ]  // bam is a list (of one BAM) so return just the one BAM
                multiple: meta.n_files > 1
            }
            .set { bam_to_merge }

        // Merge files if we have multiple files per sample
        SAMTOOLS_MERGE (
            bam_to_merge.multiple.map { meta, bam, _bai -> [ meta, bam ] },
            [[],[]],
            [[],[]]
        )
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

        // Combine merged with unmerged bam files
        SAMTOOLS_MERGE.out.bam
            .join(SAMTOOLS_MERGE.out.bai, failOnMismatch:true, failOnDuplicate:true)
            .concat(bam_to_merge.single)
            .map { meta, bam, bai -> [ meta - meta.subMap('n_files'), bam, bai ] }
            .set { ch_aligned_bam }

        // Publish alignments as CRAM if requested
        if (cram_output) {
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
            ch_input_bed
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
    // Hifiasm assembly and alignment to reference genome
    //
    if(!params.skip_genome_assembly) {

        //Hifiasm assembly
        ASSEMBLY(
            CONVERT_INPUT_FILES.out.fastq
        )
        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)

        ALIGN_ASSEMBLIES (
            ASSEMBLY.out.assembled_haplotypes,
            ch_fasta,
            ch_fai,
            cram_output
        )
        ch_versions = ch_versions.mix(ALIGN_ASSEMBLIES.out.versions)
    }

    //
    // Call SNVs
    //
    if(!params.skip_snv_calling) {

        //
        // Make BED intervals, to be used for parallel SNV calling
        //
        SCATTER_GENOME (
            ch_fai,
            ch_input_bed,             // BED file to scatter
            !params.target_regions,   // Make bed from fai
            !params.skip_snv_calling,
            params.snv_calling_processes
        )
        ch_versions = ch_versions.mix(SCATTER_GENOME.out.versions)

        // Combine to create a bam_bai - interval pair for each sample
        ch_bam_bai
            .combine( SCATTER_GENOME.out.bed_intervals )
            .map { meta, bam, bai, bed, intervals ->
                [ meta + [ num_intervals: intervals ], bam, bai, bed ]
            }
            .set{ ch_snv_calling_in }

        // This subworkflow calls SNVs with DeepVariant and outputs:
        // 1. A merged and normalized VCF, containing one sample with all regions, to be used in downstream subworkflows requiring SNVs.
        // 2. A merged and normalized VCF, containing one region with all samples, to be used in annotation and ranking.
        SHORT_VARIANT_CALLING (
            ch_snv_calling_in,
            ch_fasta,
            ch_fai,
            SCATTER_GENOME.out.bed,
            ch_par
        )
        ch_versions = ch_versions.mix(SHORT_VARIANT_CALLING.out.versions)

        // SNV QC
        SHORT_VARIANT_CALLING.out.snp_calls_vcf
            .join(SHORT_VARIANT_CALLING.out.snp_calls_tbi)
            .set { ch_snv_stats_in }

        BCFTOOLS_STATS ( ch_snv_stats_in, [[],[]], [[],[]], [[],[]], [[],[]], [[],[]] )
        ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.collect{it[1]}.ifEmpty([]))

        SHORT_VARIANT_CALLING.out.family_bcf
            .join( SHORT_VARIANT_CALLING.out.family_csi, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_vcf_tbi_per_region }
    }

    //
    // Annotate SNVs
    //
    if(!params.skip_snv_annotation) {

        // Annotates family VCFs per variant call region
        SNV_ANNOTATION(
            SHORT_VARIANT_CALLING.out.family_bcf,
            ch_databases.map { _meta, databases -> databases }.collect(),
            ch_fasta,
            ch_fai.map { name, fai -> [ [ id: name ], fai ] },
            PREPARE_GENOME.out.vep_resources.map { _meta, cache -> cache },
            params.vep_cache_version,
            ch_vep_plugin_files.collect(),
            (params.cadd_resources && params.cadd_prescored_indels),
            params.echtvar_snv_databases,
            ch_cadd_header,
            ch_cadd_resources,
            ch_cadd_prescored_indels
        )
        ch_versions = ch_versions.mix(SNV_ANNOTATION.out.versions)

        SNV_ANNOTATION.out.vcf
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
                params.filter_variants_hgnc_ids
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

        // If both CNV-calling and SV annotation is off, merged variants are output from here
        CALL_SVS (
            ch_bam_bai,
            params.sv_caller,
            ch_tandem_repeats
        )
        ch_versions = ch_versions.mix(CALL_SVS.out.versions)

        CALL_SVS.out.family_vcf
            .set { annotate_svs_in }
    }

    //
    // Call CNVs with HiFiCNV
    //
    if(!params.skip_cnv_calling) {

        CALL_CNVS (
            ch_bam_bai.join(SHORT_VARIANT_CALLING.out.snp_calls_vcf, failOnMismatch:true, failOnDuplicate:true),
            ch_fasta,
            ch_expected_xy_bed,
            ch_expected_xx_bed,
            ch_exclude_bed
        )
        ch_versions = ch_versions.mix(CALL_CNVS.out.versions)

    }

    //
    // Merge SVs and CNVs if we've called both SVs and CNVs
    //
    if (!params.skip_cnv_calling && !params.skip_sv_calling) {

        CALL_SVS.out.family_vcf
            .join(CALL_CNVS.out.family_vcf, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, svs, cnvs -> [ meta, [ svs, cnvs ] ] }
            .set { svdb_merge_svs_cnvs_in }

        // If SV annotation is off, merged variants are output from here (should be a merge and index subworkflow)
        SVDB_MERGE_SVS_CNVS (
            svdb_merge_svs_cnvs_in,
            ['svs', 'cnvs'], // Because SVs have better breakpoint resolution, give them priority
            true
        )
        ch_versions = ch_versions.mix(SVDB_MERGE_SVS_CNVS.out.versions)

        SVDB_MERGE_SVS_CNVS.out.vcf
            .set { annotate_svs_in }
    }

    //
    // Annotate SVs
    //
    if (!params.skip_sv_annotation) {

        // If annotation is on, then merged variants are output from here
        ANNOTATE_SVS (
            annotate_svs_in,
            ch_fasta,
            ch_svdb_sv_databases,
            PREPARE_GENOME.out.vep_resources.map { _meta, cache -> cache },
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
                params.filter_variants_hgnc_ids
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

    ch_collect_svs = params.skip_sv_annotation ? annotate_svs_in :
        params.skip_rank_variants ? ANN_CSQ_PLI_SVS.out.vcf :
        RANK_VARIANTS_SVS.out.vcf

    // Publish all SVs from here
    BCFTOOLS_VIEW (
        ch_collect_svs.map { meta, vcf -> [ meta, vcf, [] ] },
        [],
        [],
        []
    )

    //
    // Phase SNVs and INDELs
    //
    if(!params.skip_phasing) {

        PHASING (
            SHORT_VARIANT_CALLING.out.snp_calls_vcf,
            SHORT_VARIANT_CALLING.out.snp_calls_tbi,
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
            ch_input_bed,
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
