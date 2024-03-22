/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME             } from '../subworkflows/local/prepare_genome'
include { ASSEMBLY                   } from '../subworkflows/local/genome_assembly'
include { ASSEMBLY_VARIANT_CALLING   } from '../subworkflows/local/assembly_variant_calling'
include { ALIGN_READS                } from '../subworkflows/local/align_reads'
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

// nf-core
include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_skierfe_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SKIERFE {

    take:
    ch_sample

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Mandatory input files
    ch_fasta          = Channel.fromPath(params.fasta).map { it -> [it.simpleName, it] }.collect()

    // Optional input files
    ch_extra_snfs      = params.extra_snfs      ? Channel.fromSamplesheet('extra_snfs' , immutable_meta: false)
                                                : Channel.empty()
    ch_extra_gvcfs     = params.extra_gvcfs     ? Channel.fromSamplesheet('extra_gvcfs', immutable_meta: false)
                                                : Channel.empty()
    ch_tandem_repeats  = params.tandem_repeats  ? Channel.fromPath(params.tandem_repeats).collect()
                                                : Channel.value([])
    ch_bed             = params.bed             ? Channel.fromPath(params.bed).map{ [ it.getSimpleName(), it]}.collect()
                                                : Channel.empty()
    ch_input_bed       = params.bed             ? Channel.fromPath(params.bed).map{ [ it.getSimpleName(), it]}.collect()
                                                : Channel.value([])

    // Conditional input files that has to be set depending on which workflow is run
    ch_par             = params.dipcall_par     ? Channel.fromPath(params.dipcall_par).collect()
                                                : ''
    ch_trgt_bed        = params.trgt_repeats    ? Channel.fromPath(params.trgt_repeats).collect()
                                                : ''
    ch_databases       = params.snp_db          ? Channel.fromSamplesheet('snp_db', immutable_meta: false).map{it[1]}.collect()
                                                : ''
    ch_vep_cache       = params.ch_vep_cache    ? Channel.fromPath(params.vep_cache).collect()
                                                : ''
    ch_expected_xy_bed = params.hificnv_xy      ? Channel.fromPath(params.hificnv_xy).collect()
                                                : ''
    ch_expected_xx_bed = params.hificnv_xx      ? Channel.fromPath(params.hificnv_xx).collect()
                                                : ''
    ch_exclude_bed     = params.hificnv_exclude ? Channel.fromPath(params.hificnv_exclude).collect()
                                                : ''

    if(!params.skip_qc) {

        // Fastq QC
        FASTQC( ch_sample )
        FQCRS( ch_sample )

        // Gather versions
        ch_versions = ch_versions.mix(FASTQC.out.versions)
        ch_versions = ch_versions.mix(FQCRS.out.versions)
    }

    // Index genome
    PREPARE_GENOME( ch_fasta )
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

        BUILD_INTERVALS.out.bed
            .set{ ch_bed }
    }

    // Assembly workflow
    if(!params.skip_assembly_wf) {

        //Hifiasm assembly
        ASSEMBLY( ch_sample )
        // Run dipcall
        ASSEMBLY_VARIANT_CALLING( ASSEMBLY.out.assembled_haplotypes, fasta, fai , ch_par)

        // Gather versions
        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
        ch_versions = ch_versions.mix(ASSEMBLY_VARIANT_CALLING.out.versions)
    }

    if(!params.skip_mapping_wf) {

        ALIGN_READS( ch_sample, mmi)

        bam     = ALIGN_READS.out.bam
        bai     = ALIGN_READS.out.bai
        bam_bai = ALIGN_READS.out.bam_bai


        // TODO: parallel_snv should only be allowed when snv calling is active
        // TODO: move inside PREPARE GENOME, but only run if(parallel_snv > 1)
        // Split BED/Genome into equal chunks
        // 13 is a good number since no bin is larger than chr1 & it will not overload SLURM

        SPLIT_BED_CHUNKS(ch_bed, params.parallel_snv)

        // Combine to create a bam_bai - chunk pair for each sample
        // Do this here, pre-process or inside SNV-calling?
        bam_bai
            .combine(SPLIT_BED_CHUNKS.out
                    .split_beds
                    .flatten())
            .set{ ch_snv_calling_in }

        QC_ALIGNED_READS( bam_bai, fasta, ch_input_bed )

        // Call SVs with Sniffles2
        STRUCTURAL_VARIANT_CALLING( bam_bai , ch_extra_snfs, fasta, fai, ch_tandem_repeats )

        // Gather versions
        ch_versions = ch_versions.mix(ALIGN_READS.out.versions)
        ch_versions = ch_versions.mix(QC_ALIGNED_READS.out.versions)
        ch_versions = ch_versions.mix(STRUCTURAL_VARIANT_CALLING.out.versions)

        if(!params.skip_short_variant_calling) {
            // Call SNVs with DeepVariant/DeepTrio
            SHORT_VARIANT_CALLING( ch_snv_calling_in , ch_extra_gvcfs, fasta, fai, ch_bed )
            ch_versions = ch_versions.mix(SHORT_VARIANT_CALLING.out.versions)

            if(!params.skip_snv_annotation) {
                SNV_ANNOTATION(SHORT_VARIANT_CALLING.out.combined_bcf, SHORT_VARIANT_CALLING.out.snp_calls_vcf, ch_databases, fasta, ch_vep_cache)
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
                hap_bam_bai = PHASING.out.haplotagged_bam_bai

                // Gather versions
                ch_versions = ch_versions.mix(PHASING.out.versions)


                if(!params.skip_methylation_wf) {
                    // Pileup methylation with modkit
                    METHYLATION( hap_bam_bai, fasta, fai, ch_bed )

                    // Gather versions
                    ch_versions = ch_versions.mix(METHYLATION.out.versions)
                }

                if(!params.skip_repeat_wf) {
                    // Repeat analysis with TRGT

                    // Hack read names
                    if (params.preset == "ONT_R10") {
                        CONVERT_ONT_READ_NAMES(hap_bam_bai)
                        ch_repeat_analysis_in = CONVERT_ONT_READ_NAMES.out.bam_bai
                    } else {
                        ch_repeat_analysis_in = hap_bam_bai
                    }

                    REPEAT_ANALYSIS( ch_repeat_analysis_in, fasta, fai, ch_trgt_bed )

                    // Gather versions
                    ch_versions = ch_versions.mix(REPEAT_ANALYSIS.out.versions)
                }
            }
        }
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //

    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))
    ch_multiqc_files                      = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

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
