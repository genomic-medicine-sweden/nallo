/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//validateParameters()

WorkflowSkierfe.initialise(params, log)

def checkPathParamList = [ params.input, 
                           params.multiqc_config, 
                           params.fasta, 
                           params.ped, 
                           params.extra_snfs, 
                           params.extra_gvcfs,
                           params.dipcall_par,
                           params.tandem_repeats,
                           params.trgt_repeats
                         ]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Chech mandatory input files
ch_sample         = Channel.fromSamplesheet('input', immutable_meta: true)
ch_fasta          = Channel.fromPath(params.fasta).map { it -> [it.simpleName, it] }.collect()

// Check optional input files
ch_extra_snfs     = params.extra_snfs     ? Channel.fromSamplesheet('extra_snfs' , immutable_meta: false) : Channel.empty()
ch_extra_gvcfs    = params.extra_gvcfs    ? Channel.fromSamplesheet('extra_gvcfs', immutable_meta: false) : Channel.empty()

ch_tandem_repeats = params.tandem_repeats ? Channel.fromPath(params.tandem_repeats).collect()             : Channel.empty()

// Check mandatory files if not skipping files
if (!params.skip_assembly_wf) { 
    if (params.ped) { ch_input_ped = Channel.fromPath(params.ped) } else { exit 1, 'Input PED-file not specified - needed for dipcall!' }
    if (params.dipcall_par) { ch_par = Channel.fromPath(params.dipcall_par).collect() } else { exit 1, 'Input PAR-file not specified!' }
}

if(params.variant_caller == 'deeptrio' | params.hifiasm_mode == 'trio-binning' | !params.skip_assembly_wf | !params.skip_repeat_wf) {
    if (params.ped) { ch_input_ped = Channel.fromPath(params.ped) } else { exit 1, 'Input PED-file not specified - needed for trios, dipcall and repeats!' }
}

if (!params.skip_repeat_wf) {
    if(!params.skip_short_variant_calling) {
        if (params.trgt_repeats) { ch_trgt_bed = Channel.fromPath(params.trgt_repeats).collect() } else { exit 1, 'TGT repeat BED not specified.'}
    } else {
        exit 1, 'Short variant calling required for repeat analysis'
    }
}

// Validate workflows for different presets
validPresets = ["revio", "pacbio", "ONT_R9", "ONT_R10"]

def getValidCallers(preset) {
    switch(preset) {
        case "revio":
            return ["deepvariant"]
        case "pacbio":
            return ["deepvariant", "deeptrio", "pepper_margin_deepvariant"]
        case "ONT_R9":
            return ["pepper_margin_deepvariant"]
        case "ONT_R10":
            return ["deepvariant", "pepper_margin_deepvariant"]
    }
}

def getValidWorkflows(preset) {
    switch(preset) {
        case "pacbio":
            return ["skip_methylation_wf"]
        case "ONT_R9":
            return ["skip_repeat_wf"]
        case "ONT_R10":
            return ["skip_repeat_wf"]
    } 
}

if(params.preset !in validPresets) {
    exit 1, "Valid presets are: $validPresets"
} else if(params.preset == "pacbio" & !params.skip_methylation_wf) { 
    exit 1, "Preset \'$params.preset\' cannot be run without: " + getValidWorkflows(params.preset)
} else if(params.preset == "ONT_R9" & !params.skip_repeat_wf) { 
    exit 1, "Preset \'$params.preset\' cannot be run without: " + getValidWorkflows(params.preset)
} else if(params.preset == "ONT_R10" & !params.skip_repeat_wf) { 
    exit 1, "Preset \'$params.preset\' cannot be run without: " + getValidWorkflows(params.preset)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { FQCRS                            } from '../modules/local/fqcrs'

include { PED_CHECK                        } from '../subworkflows/local/ped_check'
include { PREPARE_GENOME                   } from '../subworkflows/local/prepare_genome'
include { ASSEMBLY                         } from '../subworkflows/local/genome_assembly'
include { ASSEMBLY_VARIANT_CALLING         } from '../subworkflows/local/assembly_variant_calling'
include { ALIGN_READS                      } from '../subworkflows/local/align_reads'
include { QC_ALIGNED_READS                 } from '../subworkflows/local/qc_aligned_reads'
include { STRUCTURAL_VARIANT_CALLING       } from '../subworkflows/local/structural_variant_calling'
include { SHORT_VARIANT_CALLING            } from '../subworkflows/local/short_variant_calling'
include { REPEAT_ANALYSIS                  } from '../subworkflows/local/repeat_analysis'
include { METHYLATION                      } from '../subworkflows/local/methylation'
include { PHASING                          } from '../subworkflows/local/phasing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MOSDEPTH                    } from '../modules/nf-core/mosdepth/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SKIERFE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet(s), validate and stage input files
    //
    
    if(params.variant_caller == 'deeptrio' | params.hifiasm_mode == 'trio-binning' | !params.skip_assembly_wf | !params.skip_repeat_wf) { 
      PED_CHECK(ch_input_ped)
        .ch_ped_processed
        .set{ ch_ped }
    } else {
        ch_ped = Channel.empty()
    }
   
    // Fastq QC 
    FASTQC( ch_sample )
    FQCRS ( ch_sample )
    
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_versions = ch_versions.mix(FQCRS.out.versions)
    
    // Index genome 
    PREPARE_GENOME( ch_fasta )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    
    // Gather indices
    fasta = ch_fasta
    fai   = PREPARE_GENOME.out.fai
    mmi   = PREPARE_GENOME.out.mmi
    
    // Assembly workflow 
    if(!params.skip_assembly_wf) {

        //Hifiasm assembly
        ASSEMBLY( ch_sample, ch_ped )
        // Run dipcall
        ASSEMBLY_VARIANT_CALLING( ASSEMBLY.out.assembled_haplotypes, fasta, fai, ch_ped, ch_par )
        
        // Gather versions 
        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
        ch_versions = ch_versions.mix(ASSEMBLY_VARIANT_CALLING.out.versions)
    }

    if(!params.skip_mapping_wf) {
        
        // Align reads
        ALIGN_READS( ch_sample, mmi )
        
        bam_bai = ALIGN_READS.out.bam_bai
        bam     = ALIGN_READS.out.bam
        bai     = ALIGN_READS.out.bai
        
        QC_ALIGNED_READS( bam, bai, fasta )
         
        // Call SVs with Sniffles2
        if(params.tandem_repeats) {
            STRUCTURAL_VARIANT_CALLING( bam_bai , ch_extra_snfs, fasta, fai, ch_tandem_repeats )
        } else {
            STRUCTURAL_VARIANT_CALLING( bam_bai , ch_extra_snfs, fasta, fai, [] )
        }
        
        // Gather versions 
        ch_versions = ch_versions.mix(STRUCTURAL_VARIANT_CALLING.out.versions)
        ch_versions = ch_versions.mix(ALIGN_READS.out.versions)
        ch_versions = ch_versions.mix(QC_ALIGNED_READS.out.versions)
        
        if(!params.skip_short_variant_calling) {
            // Call SNVs with DeepVariant/DeepTrio
            SHORT_VARIANT_CALLING( bam_bai, ch_extra_gvcfs, fasta, fai, ch_ped )
            ch_versions = ch_versions.mix(SHORT_VARIANT_CALLING.out.versions)
            // TODO: if !params.skip_structural_variant_calling
            PHASING ( SHORT_VARIANT_CALLING.out.snp_calls_vcf, STRUCTURAL_VARIANT_CALLING.out.ch_sv_calls_vcf, bam_bai, fasta, fai)
            ch_versions = ch_versions.mix(PHASING.out.versions)
            
            if(!params.skip_methylation_wf) {
                METHYLATION( WHATSHAP.out.haplotagged_bam_bai, SHORT_VARIANT_CALLING.out.snp_calls_vcf, fasta, fai )
                ch_versions = ch_versions.mix(METHYLATION.out.versions)
            }
        
            if(!params.skip_repeat_wf) {
                REPEAT_ANALYSIS( PHASING.out.haplotagged_bam_bai, fasta, fai, ch_ped, ch_trgt_bed)
            }
        }
        
       

    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //

    workflow_summary    = WorkflowSkierfe.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSkierfe.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
