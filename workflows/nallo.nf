include { fromSamplesheet } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME             } from '../subworkflows/local/prepare_genome'
include { BAM_TO_FASTQ               } from '../subworkflows/local/bam_to_fastq'
include { BAM_INFER_SEX              } from '../subworkflows/local/bam_infer_sex'
include { ASSEMBLY                   } from '../subworkflows/local/genome_assembly'
include { ASSEMBLY_VARIANT_CALLING   } from '../subworkflows/local/assembly_variant_calling'
include { QC_ALIGNED_READS           } from '../subworkflows/local/qc_aligned_reads'
include { STRUCTURAL_VARIANT_CALLING } from '../subworkflows/local/structural_variant_calling'
include { SHORT_VARIANT_CALLING      } from '../subworkflows/local/short_variant_calling'
include { CNV                        } from '../subworkflows/local/cnv'
include { REPEAT_ANALYSIS            } from '../subworkflows/local/repeat_analysis'
include { METHYLATION                } from '../subworkflows/local/methylation'
include { PHASING                    } from '../subworkflows/local/phasing'
include { SNV_ANNOTATION             } from '../subworkflows/local/snv_annotation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// local
include { FQCRS                  } from '../modules/local/fqcrs'
include { CONVERT_ONT_READ_NAMES } from '../modules/local/convert_ont_read_names'
include { BUILD_INTERVALS        } from '../modules/local/build_intervals/main'
include { SPLIT_BED_CHUNKS       } from '../modules/local/split_bed_chunks/main'
include { SAMTOOLS_MERGE         } from '../modules/nf-core/samtools/merge/main'

// nf-core
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nallo_pipeline'

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

    // Mandatory input files
    ch_fasta          = Channel.fromPath(params.fasta).map { it -> [it.simpleName, it] }.collect()

    // Optional input files
    ch_extra_snfs      = params.extra_snfs      ? Channel.fromSamplesheet('extra_snfs')
                                                : Channel.empty()
    ch_extra_gvcfs     = params.extra_gvcfs     ? Channel.fromSamplesheet('extra_gvcfs')
                                                : Channel.empty()
    ch_tandem_repeats  = params.tandem_repeats  ? Channel.fromPath(params.tandem_repeats).map{ [ it.getSimpleName(), it]}.collect()
                                                : Channel.value([[],[]])
    ch_bed             = params.bed             ? Channel.fromPath(params.bed).map{ [ it.getSimpleName(), it]}.collect()
                                                : Channel.empty()
    ch_input_bed       = params.bed             ? Channel.fromPath(params.bed).map{ [ it.getSimpleName(), it]}.collect()
                                                : Channel.value([[],[]])

    // Conditional input files that has to be set depending on which workflow is run
    ch_par             = params.dipcall_par     ? Channel.fromPath(params.dipcall_par).collect()
                                                : ''
    ch_trgt_bed        = params.trgt_repeats    ? Channel.fromPath(params.trgt_repeats).collect()
                                                : ''
    ch_databases       = params.snp_db          ? Channel.fromSamplesheet('snp_db', immutable_meta: false).map{it[1]}.collect()
                                                : ''
    ch_vep_cache_unprocessed = params.vep_cache ? Channel.fromPath(params.vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                : ''
    ch_expected_xy_bed = params.hificnv_xy      ? Channel.fromPath(params.hificnv_xy).collect()
                                                : ''
    ch_expected_xx_bed = params.hificnv_xx      ? Channel.fromPath(params.hificnv_xx).collect()
                                                : ''
    ch_exclude_bed     = params.hificnv_exclude ? Channel.fromPath(params.hificnv_exclude).collect()
                                                : ''
    ch_somalier_sites  = params.somalier_sites  ? Channel.fromPath(params.somalier_sites).map { [it.getSimpleName(), it ] }.collect()
                                                : Channel.value([[],[]])

    // Check parameter that doesn't conform to schema validation here
    if (params.split_fastq != 0 && (params.split_fastq < 2 || params.split_fastq > 999 )) { exit 1, '--split_fastq must be 0, or between 2 and 999'}
    if (params.parallel_snv == 0 ) { exit 1, '--parallel_snv must be > 0'}

    // Create PED from samplesheet
    ch_pedfile = ch_input.toList().map { file(CustomFunctions.makePed(it, params.outdir)) }

    //
    // Main workflow
    //
    BAM_TO_FASTQ ( ch_input )
    ch_versions = ch_versions.mix(BAM_TO_FASTQ.out.versions)

    BAM_TO_FASTQ.out.fastq
        .set { ch_sample }

    // Now this will be done per file and not per sample, not ideal
    if(!params.skip_qc) {

        // Fastq QC
        FASTQC( ch_sample )
        ch_versions = ch_versions.mix(FASTQC.out.versions)

        FQCRS( ch_sample )
        ch_versions = ch_versions.mix(FQCRS.out.versions)
    }

    // Index genome
    PREPARE_GENOME( ch_fasta, ch_vep_cache_unprocessed )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Gather indices
    fasta = PREPARE_GENOME.out.fasta
    fai   = PREPARE_GENOME.out.fai
    mmi   = PREPARE_GENOME.out.mmi

    // Move this inside prepare genome?

    // If no BED-file is provided then build intervals from reference
    if(!params.bed) {
        fai
            .map{ name, fai -> [['id':name], fai] }
            .set{ ch_build_intervals_in }

        BUILD_INTERVALS( ch_build_intervals_in )
        ch_versions = ch_versions.mix(BUILD_INTERVALS.out.versions)

        BUILD_INTERVALS.out.bed
            .set{ ch_bed }
    }

    if(!params.skip_mapping_wf) {

        // Split fastq
        if (params.split_fastq > 0) {

            FASTP( ch_sample, [], [], [] )
            ch_versions = ch_versions.mix(FASTP.out.versions)

            reads_for_alignment = FASTP.out.reads.transpose()

        } else {
            reads_for_alignment = ch_sample
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
            .groupTuple()
            .branch { meta, bam, bai ->
                single:   meta.n_files <= 1
                    return [ meta, bam[0], bai[0] ]  // bam is a list (of one BAM) so return just the one BAM
                multiple: meta.n_files > 1
            }
            .set { bam_to_merge }

        // Merge files if we have mutiple files per sample
        SAMTOOLS_MERGE( bam_to_merge.multiple.map { meta, bam, bai -> [ meta, bam ] }, [[],[]], [[],[]], 'bai' )
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

        // Combine merged with unmerged bams
        SAMTOOLS_MERGE.out.bam
            .join(SAMTOOLS_MERGE.out.index)
            .concat( bam_to_merge.single )
            .set { bam_infer_sex_in }

        // Infer sex if sex unknown
        BAM_INFER_SEX ( bam_infer_sex_in, fasta, fai, ch_somalier_sites, ch_pedfile )
        ch_versions = ch_versions.mix(BAM_INFER_SEX.out.versions)

        bam     = BAM_INFER_SEX.out.bam
        bai     = BAM_INFER_SEX.out.bai
        bam_bai = BAM_INFER_SEX.out.bam_bai

        // Assembly workflow
        if(!params.skip_assembly_wf) {

            //Hifiasm assembly
            ASSEMBLY( ch_sample )
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

        // TODO: parallel_snv should only be allowed when snv calling is active
        // TODO: move inside PREPARE GENOME, but only run if(parallel_snv > 1)
        // Split BED/Genome into equal chunks
        // 13 is a good number since no bin is larger than chr1 & it will not overload SLURM

        SPLIT_BED_CHUNKS(ch_bed, params.parallel_snv)
        ch_versions = ch_versions.mix(SPLIT_BED_CHUNKS.out.versions)

        // Combine to create a bam_bai - chunk pair for each sample
        // Do this here, pre-process or inside SNV-calling?
        bam_bai
            .combine(SPLIT_BED_CHUNKS.out
                    .split_beds
                    .flatten())
            .set{ ch_snv_calling_in }

        QC_ALIGNED_READS( bam_bai, fasta, ch_input_bed )
        ch_versions = ch_versions.mix(QC_ALIGNED_READS.out.versions)

        // Call SVs with Sniffles2
        STRUCTURAL_VARIANT_CALLING( bam_bai , ch_extra_snfs, fasta, fai, ch_tandem_repeats )
        ch_versions = ch_versions.mix(STRUCTURAL_VARIANT_CALLING.out.versions)

        if(!params.skip_short_variant_calling) {
            // Call SNVs with DeepVariant/DeepTrio
            SHORT_VARIANT_CALLING( ch_snv_calling_in , ch_extra_gvcfs, fasta, fai, ch_bed )
            ch_versions = ch_versions.mix(SHORT_VARIANT_CALLING.out.versions)

            if(!params.skip_snv_annotation) {

                def ch_vep_cache

                if (params.vep_cache) {
                    if (params.vep_cache.endsWith("tar.gz")) {
                        ch_vep_cache = PREPARE_GENOME.out.vep_resources
                    } else {
                        ch_vep_cache = Channel.fromPath(params.vep_cache).collect()
                    }
                } else {
                        ch_vep_cache = Channel.value([])
                }

                SNV_ANNOTATION(
                    SHORT_VARIANT_CALLING.out.combined_bcf,
                    SHORT_VARIANT_CALLING.out.snp_calls_vcf,
                    ch_databases,
                    fasta,
                    ch_vep_cache,
                    params.vep_cache_version
                )
                ch_versions = ch_versions.mix(SNV_ANNOTATION.out.versions)
            }

            if(params.preset != 'ONT_R10') {

                bam_bai
                    .join(SHORT_VARIANT_CALLING.out.snp_calls_vcf)
                    .groupTuple()
                    .set { cnv_workflow_in }

                if(!params.skip_cnv_calling) {
                    CNV(cnv_workflow_in, fasta, ch_expected_xy_bed, ch_expected_xx_bed, ch_exclude_bed)
                    ch_versions = ch_versions.mix(CNV.out.versions)
                }
            }



            if(!params.skip_phasing_wf) {
                // Phase variants with WhatsHap
                PHASING( SHORT_VARIANT_CALLING.out.snp_calls_vcf, STRUCTURAL_VARIANT_CALLING.out.ch_sv_calls_vcf, bam_bai, fasta, fai)
                ch_versions = ch_versions.mix(PHASING.out.versions)

                hap_bam_bai = PHASING.out.haplotagged_bam_bai

                if(!params.skip_methylation_wf) {
                    // Pileup methylation with modkit
                    METHYLATION( hap_bam_bai, fasta, fai, ch_bed )
                    ch_versions = ch_versions.mix(METHYLATION.out.versions)
                }

                if(!params.skip_repeat_wf) {
                    // Repeat analysis with TRGT

                    // Hack read names
                    if (params.preset == "ONT_R10") {
                        CONVERT_ONT_READ_NAMES(hap_bam_bai)
                        ch_versions = ch_versions.mix(CONVERT_ONT_READ_NAMES.out.versions)

                        ch_repeat_analysis_in = CONVERT_ONT_READ_NAMES.out.bam_bai
                    } else {
                        ch_repeat_analysis_in = hap_bam_bai
                    }

                    REPEAT_ANALYSIS( ch_repeat_analysis_in, fasta, fai, ch_trgt_bed )
                    ch_versions = ch_versions.mix(REPEAT_ANALYSIS.out.versions)
                }
            }
        }
    }

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_INFER_SEX.out.somalier_samples.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_INFER_SEX.out.somalier_pairs.map{it[1]}.collect().ifEmpty([]))

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
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

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
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
