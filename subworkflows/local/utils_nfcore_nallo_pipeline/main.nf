//
// Subworkflow with functionality specific to the genomic-medicine-sweden/nallo pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE DEPENDENCIES (FILES AND WORKFLOWS) FOR OTHER WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// nf-validation does not support contitional file and params validation,
// add these here.
//

//
// Define subworkflows and their associated "--skip"
//
def workflowSkips = [
    assembly         : "skip_assembly_wf",
    raw_read_qc      : "skip_raw_read_qc",
    aligned_read_qc  : "skip_aligned_read_qc",
    mapping          : "skip_mapping_wf",
    snv_calling      : "skip_short_variant_calling",
    snv_annotation   : "skip_snv_annotation",
    call_paralogs    : "skip_call_paralogs",
    cnv_calling      : "skip_cnv_calling",
    phasing          : "skip_phasing_wf",
    rank_variants    : "skip_rank_variants",
    repeat_calling   : "skip_repeat_calling",
    repeat_annotation: "skip_repeat_annotation",
    methylation      : "skip_methylation_wf",
]

//
//  E.g., the CNV-calling workflow depends on mapping and snv_calling and can't run without them.
//
def workflowDependencies = [
    alinged_read_qc  : ["mapping"],
    assembly         : ["mapping"],
    call_paralogs    : ["mapping"],
    snv_calling      : ["mapping"],
    snv_annotation   : ["mapping", "snv_calling"],
    cnv_calling      : ["mapping", "snv_calling"],
    phasing          : ["mapping", "snv_calling"],
    rank_variants    : ["mapping", "snv_calling", "snv_annotation"],
    repeat_calling   : ["mapping", "snv_calling", "phasing"],
    repeat_annotation: ["mapping", "snv_calling", "phasing", "repeat_calling"],
    methylation      : ["mapping", "snv_calling", "phasing"],
]

//
// E.g., the dipcall_par file is required by the assembly workflow and the assembly workflow can't run without dipcall_par
//
def fileDependencies = [
    mapping          : ["fasta", "somalier_sites"],
    assembly         : ["fasta", "dipcall_par"], // The assembly workflow should be split into two - assembly and variant calling (requires ref)
    snv_annotation   : ["snp_db", "vep_cache", "reduced_penetrance", "score_config_snv", "variant_consequences_snv"],
    cnv_calling      : ["hificnv_xy", "hificnv_xx", "hificnv_exclude"],
    repeat_calling   : ["trgt_repeats"],
    repeat_annotation: ["variant_catalog"],
]

def parameterStatus = [
    workflow: [
        skip_short_variant_calling: params.skip_short_variant_calling,
        skip_phasing_wf           : params.skip_phasing_wf,
        skip_methylation_wf       : params.skip_methylation_wf,
        skip_rank_variants        : params.skip_rank_variants,
        skip_repeat_calling       : params.skip_repeat_calling,
        skip_repeat_annotation    : params.skip_repeat_annotation,
        skip_snv_annotation       : params.skip_snv_annotation,
        skip_call_paralogs        : params.skip_call_paralogs,
        skip_cnv_calling          : params.skip_cnv_calling,
        skip_mapping_wf           : params.skip_mapping_wf,
        skip_aligned_read_qc      : params.skip_aligned_read_qc,
        skip_raw_read_qc          : params.skip_raw_read_qc,
        skip_assembly_wf          : params.skip_assembly_wf,
    ],
    files: [
        dipcall_par             : params.dipcall_par,
        snp_db                  : params.snp_db,
        somalier_sites          : params.somalier_sites,
        vep_cache               : params.vep_cache,
        hificnv_xy              : params.hificnv_xy,
        hificnv_xx              : params.hificnv_xx,
        hificnv_exclude         : params.hificnv_exclude,
        fasta                   : params.fasta,
        trgt_repeats            : params.trgt_repeats,
        variant_catalog         : params.variant_catalog,
        score_config_snv        : params.score_config_snv,
        reduced_penetrance      : params.reduced_penetrance,
        score_config_snv        : params.score_config_snv,
        variant_consequences_snv: params.variant_consequences_snv,
    ]
]

/*
    SUBWORKFLOW TO INITIALISE PIPELINE
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters(parameterStatus, workflowSkips, workflowDependencies, fileDependencies)

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .map { meta, reads ->
            [ meta.id, meta, reads ] // add sample as groupTuple key
        }
        .map {
            validateInputSamplesheet(it)
        }
        .groupTuple() // group by sample
        .map { sample, metas, reads ->
            // Add number of files per sample _after_ splitting to meta
            [ sample, metas[0] + [n_files: metas.size() + metas.size() * Math.max(0, params.split_fastq - 1), single_end:true ], reads ]
        }
        // Convert back to [ meta, reads ]
        .flatMap {
            sample, meta, reads ->
                reads.collect { return [ meta, it ] }
        }
        .set { ch_samplesheet }

        // Check that there's samples with affected phenotype if we are ranking variants
        ch_samplesheet
            .filter { meta, reads -> meta.phenotype == 2 }
            .ifEmpty {
                if(!params.skip_rank_variants) {
                    error("No samples in samplesheet has affected phenotype (=2), --skip_rank_variants has to be active.")
                }
            }


    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//

def validateInputParameters(statusMap, workflowMap, workflowDependencies, fileDependencies) {
    genomeExistsError()
    validateParameterCombinations(statusMap, workflowMap, workflowDependencies, fileDependencies)
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {

    return input
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {

    def repeat_annotation_text = []
    def preprocessing_text     = []
    def other_citation_text    = []

    if (!params.skip_repeat_annotation) {
        repeat_annotation_text = [
            "stranger (Nilsson & Magnusson, 2021),"
        ]
    }
    preprocessing_text = [
        "FastQC (Andrews 2010),",
    ]
    other_citation_text = [
        "MultiQC (Ewels et al. 2016),",
        "."
    ]
    def concat_text = repeat_annotation_text +
        preprocessing_text +
        other_citation_text

    def citation_text = [ "Tools used in the workflow included:" ] + concat_text.unique(false) { a, b -> a <=> b } - ""
    return citation_text.join(' ').trim()
}

def toolBibliographyText() {
    def repeat_annotation_text = []
    def preprocessing_text     = []
    def other_citation_text    = []

    if (!params.skip_repeat_annotation) {
        repeat_annotation_text = [
            "<li>Nilsson, D., & Magnusson, M. (2021). Moonso/stranger v0.9.1 (v0.9.1) [Computer software]. Zenodo. https://zenodo.org/doi/10.5281/zenodo.3841097</li>"
        ]
    }

    preprocessing_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
    ]

    other_citation_text = [
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
    ].join(' ').trim()

    def concat_text = repeat_annotation_text +
        preprocessing_text +
        other_citation_text

    def reference_text = concat_text.unique(false) { a, b -> a <=> b } - ""
    return reference_text.join(' ').trim()
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Validate  workflow skip combinations
//
def validateParameterCombinations(statusMap, workflowMap, workflowDependencies, fileDependencies) {
    // Array to store errors
    def errors = []
    // For each of the "workflow", "files"
    statusMap.each { paramsType, allParams ->
        // Go through all params and their status
        statusMap[paramsType].each { param, paramStatus ->
            switch (paramsType) {
                case "files":
                    checkFileDependencies(param, fileDependencies, statusMap, workflowMap, errors)
                    break
                case "workflow":
                    checkWorkflowDependencies(param, workflowDependencies, statusMap, workflowMap, errors)
                    break
                default:
                    break
            }
        }
    }
    // Give error if there are any
    if(errors) {
        def error_string =
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  " + errors.join("\n  ") + "\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Lookup all workflows that needs to be active for another workflow
//
def checkWorkflowDependencies(String skip, Map combinationsMap, Map statusMap, Map workflowMap, List errors) {

    // Lookup the workflow associated with the --skip_xxx parameter
    currentWorkflow = workflowMap.find { key, mapValue -> mapValue == skip }?.key

    // If the --skip is not set, then the workflow is active, give no error
    workflowIsActive = !statusMap["workflow"][skip]
    if(workflowIsActive) {
        return
    }

    // Get all other worflows that are required for a certain workflow
    def requiredWorkflows = combinationsMap.findAll { it.value.contains(currentWorkflow) }.keySet()
    // If no direct dependencies are found or combinationsMap does not contain the workflow, return an empty list
    if (!requiredWorkflows) {
        return []
    }
    // Collect the required --skips that are not active for the current workflow
    def dependencyString = findRequiredSkips("workflow", requiredWorkflows, statusMap, workflowMap)
        .collect { [ '--', it ].join('') }
        .join(" ")
    // If all reqired sets are set, give no error
    if (!dependencyString) {
        return
    }
    errors << "--$skip is active, the pipeline has to be run with: $dependencyString"
    return errors
}

//
// Lookup if a file is required by any workflows, and add to errors
//
def checkFileDependencies(String file, Map combinationsMap, Map statusMap, Map workflowMap, List errors) {
    // Get the the workflow required by file
    def workflowThatRequiresFile = findKeyForValue(file, combinationsMap)
    // Get the "--skip" for that workflow
    def workflowSkip = workflowMap[workflowThatRequiresFile]
    // Get the status of the "--skip", if false then workflow is active
    def WorkflowIsActive = !statusMap["workflow"][workflowSkip]
    // Get the file path
    def FilePath = statusMap["files"][file]
    // If the workflow that requires the file is active & theres no file available
    if(WorkflowIsActive && FilePath == null) {
        errors << "--$workflowSkip is NOT active, the following files are required: --$file"
    }
    return errors
}

//
// Find the workflow skips that are not currently active
//
def findRequiredSkips(paramType, Set<String> requiredWorkflows, Map statusMap, Map workflowMap) {

    def requiredSkips = []

    for (currentWorkflow in requiredWorkflows) {
        // Get the skip associated with the workflow
        skip = workflowMap[currentWorkflow]

        workflowIsSkipped = !statusMap[paramType][skip]

        if(paramType == "workflow") {
            if(workflowIsSkipped) {
                requiredSkips << skip
            }
        }
    }
    return requiredSkips
}

def findKeyForValue(def valueToFind, Map map) {
    for (entry in map) {
        def key = entry.key
        def value = entry.value

        if (value instanceof List) {
            if (value.contains(valueToFind)) {
                return key
            }
        } else {
            if (value == valueToFind) {
                return key
            }
        }
    }
    return null // Value not found
}

