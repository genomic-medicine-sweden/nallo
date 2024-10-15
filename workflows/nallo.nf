include { samplesheetToList } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV     } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_SVS                            } from '../subworkflows/local/annotate_svs'
include { ANNOTATE_REPEAT_EXPANSIONS              } from '../subworkflows/local/annotate_repeat_expansions'
include { ASSEMBLY                                } from '../subworkflows/local/genome_assembly'
include { ASSEMBLY_VARIANT_CALLING                } from '../subworkflows/local/assembly_variant_calling'
include { CALL_SVS                                } from '../subworkflows/local/call_svs'
include { CONVERT_INPUT_FILES                     } from '../subworkflows/local/convert_input_files'
include { BAM_INFER_SEX                           } from '../subworkflows/local/bam_infer_sex'
include { CALL_PARALOGS                           } from '../subworkflows/local/call_paralogs'
include { CALL_REPEAT_EXPANSIONS                  } from '../subworkflows/local/call_repeat_expansions'
include { CNV                                     } from '../subworkflows/local/cnv'
include { METHYLATION                             } from '../subworkflows/local/methylation'
include { PHASING                                 } from '../subworkflows/local/phasing'
include { PREPARE_GENOME                          } from '../subworkflows/local/prepare_genome'
include { QC_ALIGNED_READS                        } from '../subworkflows/local/qc_aligned_reads'
include { RANK_VARIANTS as RANK_VARIANTS_SNV      } from '../subworkflows/local/rank_variants'
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
include { ECHTVAR_ENCODE                                    } from '../modules/local/echtvar/encode/main'
include { SAMTOOLS_MERGE                                    } from '../modules/nf-core/samtools/merge/main'

// nf-core
include { BCFTOOLS_CONCAT                                   } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_PLUGINSPLIT as BCFTOOLS_PLUGINSPLIT_SNVS } from '../modules/nf-core/bcftools/pluginsplit/main'
include { BCFTOOLS_PLUGINSPLIT as BCFTOOLS_PLUGINSPLIT_SVS  } from '../modules/nf-core/bcftools/pluginsplit/main'
include { BCFTOOLS_SORT                                     } from '../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_STATS                                    } from '../modules/nf-core/bcftools/stats/main'
include { MINIMAP2_ALIGN                                    } from '../modules/nf-core/minimap2/align/main'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { SPLITUBAM                                         } from '../modules/nf-core/splitubam/main'
include { paramsSummaryMap                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../subworkflows/local/utils_nfcore_nallo_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NALLO {

    take:
    ch_input

    main:
    ch_vep_cache     = Channel.value([])
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Optional input files that has to be set depending on which workflow is run
    ch_cadd_header              = Channel.fromPath("$projectDir/assets/cadd_to_vcf_header_-1.0-.txt", checkIfExists: true).collect()
    ch_cadd_resources           = params.cadd_resources             ? Channel.fromPath(params.cadd_resources).collect()
                                                                    : ''
    ch_cadd_prescored           = params.cadd_prescored             ? Channel.fromPath(params.cadd_prescored).collect()
                                                                    : ''
    ch_fasta                    = params.fasta                      ? Channel.fromPath(params.fasta).map { it -> [ it.simpleName, it ] }.collect()
                                                                    : ''
    ch_tandem_repeats           = params.tandem_repeats             ? Channel.fromPath(params.tandem_repeats).map{ [ it.simpleName, it ] }.collect()
                                                                    : Channel.value([[],[]])
    ch_input_bed                = params.bed                        ? Channel.fromPath(params.bed).map{ [ [ id:it.simpleName ] , it ] }.collect()
                                                                    : Channel.value([[],[]])
    ch_par                      = params.par_regions                ? Channel.fromPath(params.par_regions).map { [ [ id: it.simpleName ], it ] }.collect()
                                                                    : ''
    ch_trgt_bed                 = params.trgt_repeats               ? Channel.fromPath(params.trgt_repeats).map { it -> [ it.simpleName, it ] }.collect()
                                                                    : ''
    ch_variant_catalog          = params.variant_catalog            ? Channel.fromPath(params.variant_catalog).map { it -> [ it.simpleName, it ] }.collect()
                                                                    : ''
    ch_databases                = params.snp_db                     ? Channel.fromList(samplesheetToList(params.snp_db, 'assets/schema_snpdb.json')).map{ it[1] }.collect()
                                                                    : ''
    ch_variant_consequences_snv = params.variant_consequences_snv   ? Channel.fromPath(params.variant_consequences_snv).collect()
                                                                    : Channel.value([])
    ch_vep_cache_unprocessed    = params.vep_cache                  ? Channel.fromPath(params.vep_cache).map { it -> [ [ id:'vep_cache' ], it ] }.collect()
                                                                    : Channel.value([[],[]])
    ch_vep_extra_files_unsplit  = params.vep_plugin_files           ? Channel.fromPath(params.vep_plugin_files).collect()
                                                                    : ''
    ch_expected_xy_bed          = params.hificnv_xy                 ? Channel.fromPath(params.hificnv_xy).collect()
                                                                    : ''
    ch_expected_xx_bed          = params.hificnv_xx                 ? Channel.fromPath(params.hificnv_xx).collect()
                                                                    : ''
    ch_exclude_bed              = params.hificnv_exclude            ? Channel.fromPath(params.hificnv_exclude).collect()
                                                                    : ''
    ch_reduced_penetrance       = params.reduced_penetrance         ? Channel.fromPath(params.reduced_penetrance).collect()
                                                                    : Channel.value([])
    ch_score_config_snv         = params.score_config_snv           ? Channel.fromPath(params.score_config_snv).collect()
                                                                    : Channel.value([])
    ch_somalier_sites           = params.somalier_sites             ? Channel.fromPath(params.somalier_sites).map { [ it.simpleName, it ] }.collect()
                                                                    : ''
    ch_svdb_dbs                 = params.svdb_dbs                   ? Channel.fromPath(params.svdb_dbs).map { [ it.simpleName, it ] }.collect()
                                                                    : ''

    // Check parameter that doesn't conform to schema validation here
    if (params.phaser.matches('hiphase') && params.preset == 'ONT_R10') { error "The HiPhase license only permits analysis of data from PacBio. For details see: https://github.com/PacificBiosciences/HiPhase/blob/main/LICENSE.md" }

    // Read and store paths in the vep_plugin_files file
    if (params.vep_plugin_files) {
        ch_vep_extra_files_unsplit.splitCsv ( header:true )
            .map { row ->
                path = file(row.vep_files[0])
                if(path.isFile() || path.isDirectory()){
                    return [path]
                } else {
                    error("\nVep database file ${path} does not exist.")
                }
            }
            .collect()
            .set {ch_vep_extra_files}
    }

    //
    // Convert FASTQ to BAM (and vice versa if assembly workflow is active)
    //
    CONVERT_INPUT_FILES (
        ch_input,
        !params.skip_assembly_wf
    )
    ch_versions = ch_versions.mix(CONVERT_INPUT_FILES.out.versions)

    //
    // Prepare references
    //
    if(!params.skip_mapping_wf | !params.skip_assembly_wf ) {

        PREPARE_GENOME (
            ch_fasta,
            ch_vep_cache_unprocessed,
        )
        ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

        if(!params.skip_snv_annotation) {
            if (params.vep_cache) {
                if (params.vep_cache.endsWith("tar.gz")) {
                    ch_vep_cache = PREPARE_GENOME.out.vep_resources
                } else {
                    ch_vep_cache = Channel.fromPath(params.vep_cache).collect()
                }
            }
        }

        // Gather indices
        fasta = PREPARE_GENOME.out.fasta
        fai   = PREPARE_GENOME.out.fai
        mmi   = PREPARE_GENOME.out.mmi
    }

    //
    // (Split input files and), map reads to reference and merge into a single BAM per sample
    //
    if(!params.skip_mapping_wf) {

        // Split input files for alignment
        if (params.parallel_alignments > 1) {

            SPLITUBAM ( CONVERT_INPUT_FILES.out.bam )
            ch_versions = ch_versions.mix(SPLITUBAM.out.versions)

            reads_for_alignment = SPLITUBAM.out.bam.transpose()

        } else {
            reads_for_alignment = CONVERT_INPUT_FILES.out.bam
        }
        // Align (split) reads
        MINIMAP2_ALIGN ( reads_for_alignment, mmi, true, 'bai', false, false )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        // Split channel into cases where we have multiple files or single files
        MINIMAP2_ALIGN.out.bam
            .join(MINIMAP2_ALIGN.out.index)
            .map {
                meta, bam, bai ->
                    [ groupKey(meta, meta.n_files), bam, bai ]
            }
            .view()
            .groupTuple()
            .view() // Needs to be 16 files, but split_bam only produced 8 - needs error in split_ubam
            .branch { meta, bam, bai ->
                single:   meta.n_files <= 1
                    return [ meta, bam[0], bai[0] ]  // bam is a list (of one BAM) so return just the one BAM
                multiple: meta.n_files > 1
            }
            .set { bam_to_merge }

        // Merge files if we have multiple files per sample
        SAMTOOLS_MERGE ( bam_to_merge.multiple.map { meta, bam, bai -> [ meta, bam ] }, [[],[]], [[],[]], 'bai' )
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

        // Combine merged with unmerged bams
        SAMTOOLS_MERGE.out.bam
            .join(SAMTOOLS_MERGE.out.index)
            .concat(bam_to_merge.single)
            .map { meta, bam, bai -> [ meta - meta.subMap('n_files'), bam, bai ] }
            .set { bam_infer_sex_in }

        //
        // Create PED from samplesheet
        //
        ch_input
            .map { meta, files -> [ meta.project, meta ] }
            .groupTuple()
            .set { ch_samplesheet_ped_in }

        SAMPLESHEET_PED ( ch_samplesheet_ped_in )
        ch_versions = ch_versions.mix(SAMPLESHEET_PED.out.versions)

        SAMPLESHEET_PED.out.ped
            .map { project, ped -> [ [ 'id': project ], ped ] }
            .collect()
            .set { ch_samplesheet_pedfile }

        //
        // Check sex and relatedness, and update with infered sex if the sex for a sample is unknown
        //
        BAM_INFER_SEX ( bam_infer_sex_in, fasta, fai, ch_somalier_sites, ch_samplesheet_pedfile )
        ch_versions = ch_versions.mix(BAM_INFER_SEX.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(BAM_INFER_SEX.out.somalier_samples.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BAM_INFER_SEX.out.somalier_pairs.map{it[1]}.collect().ifEmpty([]))

        bam     = BAM_INFER_SEX.out.bam
        bai     = BAM_INFER_SEX.out.bai
        bam_bai = BAM_INFER_SEX.out.bam_bai

        //
        // Create PED with updated sex
        //
        bam
            .map { meta, files -> [ meta.project, meta ] }
            .groupTuple()
            .set { ch_somalier_ped_in }

        SOMALIER_PED ( ch_somalier_ped_in )
        ch_versions = ch_versions.mix(SOMALIER_PED.out.versions)

        SOMALIER_PED.out.ped
            .map { project, ped -> [ [ 'id': project ], ped ] }
            .collect()
            .set { ch_updated_pedfile }

        //
        // Run read QC with FastQC, mosdepth and cramino
        //
        if (!params.skip_qc) {

            QC_ALIGNED_READS( bam_bai, fasta, ch_input_bed )
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
            CALL_PARALOGS ( bam_bai, fasta )
            ch_versions = ch_versions.mix(CALL_PARALOGS.out.versions)
        }

        //
        // Hifiasm assembly and assembly variant calling
        //
        if(!params.skip_assembly_wf) {

            //Hifiasm assembly
            ASSEMBLY( CONVERT_INPUT_FILES.out.fastq )
            ch_versions = ch_versions.mix(ASSEMBLY.out.versions)

            // Update assembly variant calling meta with sex from somalier
            ASSEMBLY.out.assembled_haplotypes
                .map { meta, hap1, hap2 -> [ meta.id, [ hap1, hap2 ] ] }
                .set { haplotypes }

            bam
                .map { meta, bam -> [ meta.id, meta ] }
                .join( haplotypes )
                .map { id, meta, haplotypes -> [ meta, haplotypes[0], haplotypes[1] ] }
                .set { ch_assembly_variant_calling_in }

            // Run dipcall
            ASSEMBLY_VARIANT_CALLING( ch_assembly_variant_calling_in, fasta, fai , ch_par)
            ch_versions = ch_versions.mix(ASSEMBLY_VARIANT_CALLING.out.versions)
        }

        //
        // Call structural variants
        //
        CALL_SVS (
            bam_bai,
            fasta,
            fai,
            params.sv_caller,
            ch_tandem_repeats,
            ch_input_bed
        )
        ch_versions = ch_versions.mix(CALL_SVS.out.versions)

        CALL_SVS.out.ch_multisample_vcf
            .join( CALL_SVS.out.ch_multisample_tbi )
            .set { ch_split_svs_in }
        //
        // Annotate structural variants
        //
        if(!params.skip_sv_annotation) {
            ANNOTATE_SVS (
                CALL_SVS.out.ch_multisample_vcf,
                fasta,
                ch_svdb_dbs,
                ch_vep_cache,
                params.vep_cache_version,
                ch_vep_extra_files
            )

            ANNOTATE_SVS.out.vcf
                .join( ANNOTATE_SVS.out.tbi )
                .set { ch_split_svs_in }
        }

        // Split the multisample SV VCF to also publish an (annotated) VCF per sample
        BCFTOOLS_PLUGINSPLIT_SVS (
            ch_split_svs_in,
            [],
            [],
            [],
            []
        )
        ch_versions = ch_versions.mix(BCFTOOLS_PLUGINSPLIT_SVS.out.versions)

        //
        // Call (and annotate and rank) SNVs
        //
        if(!params.skip_short_variant_calling) {

            //
            // Make BED intervals, to be used for parallel SNV calling
            //
            SCATTER_GENOME (
                fai,
                ch_input_bed,
                !params.bed,
                !params.skip_short_variant_calling,
                params.parallel_snv
            )
            ch_versions = ch_versions.mix(SCATTER_GENOME.out.versions)

            // Combine to create a bam_bai - interval pair for each sample
            bam_bai
                .combine( SCATTER_GENOME.out.bed_intervals )
                .map { meta, bam, bai, bed, intervals ->
                    [ meta + [ num_intervals: intervals ], bam, bai, bed ]
                }
                .set{ ch_snv_calling_in }

            //
            // This subworkflow calls SNVs with DeepVariant and outputs:
            // 1. A merged and normalised VCF, containing one sample with all regions, to be used in downstream subworkflows requiring SNVs.
            // 2. A merged and normalised VCF, containing one region with all samples, to be used in annotation and ranking.
            //
            SHORT_VARIANT_CALLING( ch_snv_calling_in, fasta, fai, SCATTER_GENOME.out.bed, ch_par )
            ch_versions = ch_versions.mix(SHORT_VARIANT_CALLING.out.versions)

            //
            // Annotate SNVs
            //
            if(!params.skip_snv_annotation) {

                //
                // Annotates one multisample VCF per variant call region
                //
                SNV_ANNOTATION(
                    SHORT_VARIANT_CALLING.out.combined_bcf,
                    ch_databases,
                    fasta,
                    fai.map { name, fai -> [ [ id: name ], fai ] },
                    ch_vep_cache,
                    params.vep_cache_version,
                    ch_vep_extra_files,
                    (params.cadd_resources && params.cadd_prescored),
                    ch_cadd_header,
                    ch_cadd_resources,
                    ch_cadd_prescored
                )
                ch_versions = ch_versions.mix(SNV_ANNOTATION.out.versions)

                ANN_CSQ_PLI_SNV (
                    SNV_ANNOTATION.out.vcf,
                    ch_variant_consequences_snv
                )
                ch_versions = ch_versions.mix(ANN_CSQ_PLI_SNV.out.versions)

                //
                // Ranks one multisample VCF per variant call region
                //
                if(!params.skip_rank_variants) {
                    // Only run if we have affected individuals
                    RANK_VARIANTS_SNV (
                        ANN_CSQ_PLI_SNV.out.vcf_ann.filter { meta, vcf -> meta.contains_affected },
                        ch_updated_pedfile.map { meta, ped -> ped },
                        ch_reduced_penetrance,
                        ch_score_config_snv
                    )
                    ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)

                    // If there are affected individuals and RANK_VARIANTS has been run,
                    // input that to VCF concatenation
                    RANK_VARIANTS_SNV.out.vcf
                        .join( RANK_VARIANTS_SNV.out.tbi )
                        .set { ch_vcf_tbi_per_region }
                } else {
                    // otherwise grab the VCF that should have gone into RANK_VARIANTS
                    ANN_CSQ_PLI_SNV.out.vcf_ann
                        .join( ANN_CSQ_PLI_SNV.out.tbi_ann )
                        .set { ch_vcf_tbi_per_region }
                }
            } else {
                // If neither snv_annotation nor rank_variants was run, take the output from
                // SHORT_VARIANT_CALLING
                SHORT_VARIANT_CALLING.out.combined_bcf
                    .join( SHORT_VARIANT_CALLING.out.combined_csi )
                    .set { ch_vcf_tbi_per_region }
            }

            ch_vcf_tbi_per_region
                .map { meta, vcf, tbi -> [ [ id: meta.project ], vcf, tbi ] }
                .groupTuple()
                .set { ch_bcftools_concat_in }

            // Concat into a multisample VCF with all regions
            BCFTOOLS_CONCAT ( ch_bcftools_concat_in )
            ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

            // Sort and publish
            BCFTOOLS_SORT ( BCFTOOLS_CONCAT.out.vcf )
            ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

            // Make an echtvar database of all samples
            ECHTVAR_ENCODE ( BCFTOOLS_SORT.out.vcf )
            ch_versions = ch_versions.mix(ECHTVAR_ENCODE.out.versions)

            // Split multisample VCF to also publish a VCF per sample
            BCFTOOLS_PLUGINSPLIT_SNVS ( BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_SORT.out.tbi ), [], [], [], [] )
            ch_versions = ch_versions.mix(BCFTOOLS_PLUGINSPLIT_SNVS.out.versions)

            BCFTOOLS_PLUGINSPLIT_SNVS.out.vcf
                .transpose()
                .map { meta, vcf -> [ meta, vcf, [] ] }
                .set { ch_bcftools_stats_snv_in }

            BCFTOOLS_STATS ( ch_bcftools_stats_snv_in, [[],[]], [[],[]], [[],[]], [[],[]], [[],[]] )
            ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.collect{it[1]}.ifEmpty([]))

            //
            // Call CNVs with HiFiCNV
            //
            if(!params.skip_cnv_calling) {
                bam_bai
                    .join(SHORT_VARIANT_CALLING.out.snp_calls_vcf)
                    .set { cnv_workflow_in }

                CNV(cnv_workflow_in, fasta, ch_expected_xy_bed, ch_expected_xx_bed, ch_exclude_bed)
                ch_versions = ch_versions.mix(CNV.out.versions)
            }

            //
            // Phase SNVs and INDELs
            //
            if(!params.skip_phasing_wf) {

                PHASING (
                    SHORT_VARIANT_CALLING.out.snp_calls_vcf,
                    SHORT_VARIANT_CALLING.out.snp_calls_tbi,
                    bam_bai,
                    fasta,
                    fai
                )
                ch_versions = ch_versions.mix(PHASING.out.versions)

                ch_multiqc_files = ch_multiqc_files.mix(PHASING.out.stats.collect{it[1]}.ifEmpty([]))

                //
                // Create methylation pileups with modkit
                //
                if(!params.skip_methylation_wf) {

                    METHYLATION( PHASING.out.haplotagged_bam_bai, fasta, fai, ch_input_bed )
                    ch_versions = ch_versions.mix(METHYLATION.out.versions)
                }

                //
                // Call repeat expansions with TRGT
                //
                if(!params.skip_repeat_calling) {

                    CALL_REPEAT_EXPANSIONS ( PHASING.out.haplotagged_bam_bai, fasta, fai, ch_trgt_bed )
                    ch_versions = ch_versions.mix(CALL_REPEAT_EXPANSIONS.out.versions)

                    //
                    // Annotate repeat expansions with stranger
                    //
                    if(!params.skip_repeat_annotation) {
                        ANNOTATE_REPEAT_EXPANSIONS ( ch_variant_catalog, CALL_REPEAT_EXPANSIONS.out.vcf )
                        ch_versions = ch_versions.mix(ANNOTATE_REPEAT_EXPANSIONS.out.versions)
                    }
                }
            }
        }
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
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
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
