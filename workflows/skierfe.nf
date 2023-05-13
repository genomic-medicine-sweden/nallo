/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSkierfe.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.ped, params.extra_snfs, params.extra_gvcfs ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = Channel.fromPath(params.fasta) } else { exit 1, 'Input fasta not specified!' }
if (params.trio) { if (params.ped) { ch_input_ped = file(params.ped) } else { exit 1, 'Input PED-file not specified!' } }
// TODO: Should be required only if running DIPCALL
if (params.par) { ch_par = Channel.fromPath(params.par) } else { exit 1, 'Input PAR-file not specified!' }

ch_fasta = ch_fasta
  .map { it -> [it.simpleName, it] }
  .groupTuple()

// Not pretty but works for now...needs initializing?
params.extra_snfs = ''
params.extra_gvcfs = ''

// Since they are not mandatory, populate channel only if the samplesheets are provided
if (params.extra_snfs) { ch_input_snfs = file(params.extra_snfs) } else { ch_input_snfs = Channel.empty() }
if (params.extra_gvcfs) { ch_input_gvcfs = file(params.extra_gvcfs) } else { ch_input_gvcfs = Channel.empty() }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:w
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

include { INPUT_CHECK as INPUT_FASTQ_CHECK } from '../subworkflows/local/input_check'
include { INPUT_CHECK as SNFS_CHECK } from '../subworkflows/local/input_check'
include { INPUT_CHECK as GVCFS_CHECK } from '../subworkflows/local/input_check'
include { PED_CHECK } from '../subworkflows/local/ped_check'


include { ASSEMBLY } from '../subworkflows/local/genome_assembly'
include { ASSEMBLY_VARIANT_CALLING } from '../subworkflows/local/assembly_variant_calling'

include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ALIGN_READS } from '../subworkflows/local/align_reads'
include { STRUCTURAL_VARIANT_CALLING } from '../subworkflows/local/structural_variant_calling'
include { SHORT_VARIANT_CALLING } from '../subworkflows/local/short_variant_calling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
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
    ch_sample = Channel.empty()
    ch_ped = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet(s), validate and stage input files
    //

    INPUT_FASTQ_CHECK ( ch_input )
        .ch_sample.set { ch_sample }

    SNFS_CHECK ( ch_input_snfs )
        .ch_sample.set { ch_extra_snfs }
       
    GVCFS_CHECK ( ch_input_gvcfs )
        .ch_sample.set { ch_extra_gvcfs }
    
    if(params.trio) { 
      PED_CHECK(ch_input_ped)
        .ch_ped_processed
        .set{ch_ped}
    }

    ch_versions = ch_versions.mix(INPUT_FASTQ_CHECK.out.versions)

    //Hifiasm assembly
    ASSEMBLY ( ch_sample, ch_ped )

    ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
   
    //FASTQC
    FASTQC( ch_sample )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // Index the genome 
    PREPARE_GENOME ( ch_fasta )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
   
    // Run dipcall
    ASSEMBLY_VARIANT_CALLING ( ASSEMBLY.out.assembled_haplotypes, PREPARE_GENOME.out.fasta, PREPARE_GENOME.out.fai, PREPARE_GENOME.out.mmi, ch_ped, ch_par )
    ch_versions = ch_versions.mix(ASSEMBLY_VARIANT_CALLING.out.versions)

    // Align reads with pbmm2 
    ALIGN_READS ( ch_sample, PREPARE_GENOME.out.mmi )
    ch_versions = ch_versions.mix(ALIGN_READS.out.versions)
    
    // Call SVs with Sniffles2 
    STRUCTURAL_VARIANT_CALLING ( ALIGN_READS.out.bam_bai , ch_extra_snfs, ch_fasta )
    ch_versions = ch_versions.mix(STRUCTURAL_VARIANT_CALLING.out.versions)
    
    // Call SNVs with DeepVariant/DeepTrio
    SHORT_VARIANT_CALLING ( ALIGN_READS.out.bam_bai, ch_input_gvcfs, ch_fasta, PREPARE_GENOME.out.fai, ch_ped )
    ch_versions = ch_versions.mix(SHORT_VARIANT_CALLING.out.versions)

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
