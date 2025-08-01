//
// Subworkflow with functionality specific to the genomic-medicine-sweden/nallo pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
    SUBWORKFLOW TO INITIALISE PIPELINE
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
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
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // nf-validation does not support conditional file and params validation,
    // add these here.
    //

    //
    // Define subworkflows and their associated "--skip"
    //
    def workflowSkips = [
        assembly         : "skip_genome_assembly",
        mapping          : "skip_alignment",
        snv_calling      : "skip_snv_calling",
        snv_annotation   : "skip_snv_annotation",
        sv_calling       : "skip_sv_calling",
        sv_annotation    : "skip_sv_annotation",
        call_paralogs    : "skip_call_paralogs",
        peddy            : "skip_peddy",
        phasing          : "skip_phasing",
        rank_variants    : "skip_rank_variants",
        repeat_calling   : "skip_repeat_calling",
        repeat_annotation: "skip_repeat_annotation",
        methylation      : "skip_methylation_pileups",
        qc               : "skip_qc",
    ]

    //
    //  E.g., the CNV-calling workflow depends on mapping and snv_calling and can't run without them.
    //
    def workflowDependencies = [
        call_paralogs    : ["mapping"],
        snv_calling      : ["mapping"],
        qc               : ["mapping"],
        sv_calling       : ["mapping"],
        sv_annotation    : ["mapping", "sv_calling"],
        peddy            : ["mapping", "snv_calling"],
        snv_annotation   : ["mapping", "snv_calling"],
        phasing          : ["mapping", "snv_calling"],
        rank_variants    : ["mapping", "snv_calling", "snv_annotation", "sv_annotation"],
        repeat_calling   : ["mapping", "snv_calling", "phasing"],
        repeat_annotation: ["mapping", "snv_calling", "phasing", "repeat_calling"],
        methylation      : ["mapping", "snv_calling"]
    ]

    //
    // E.g., the par_regions file is required by the assembly workflow and the assembly workflow can't run without par_regions
    //
    def fileDependencies = [
        mapping          : ["fasta", "somalier_sites"],
        assembly         : ["fasta"], // The assembly workflow should perhaps be split into two - assembly and alignment (requires ref)
        snv_calling      : ["fasta", "par_regions"],
        snv_annotation   : ["vep_cache", "vep_plugin_files", "variant_consequences_snvs"],
        sv_calling       : ["fasta"],
        sv_annotation    : ["svdb_sv_databases", "vep_cache", "vep_plugin_files", "variant_consequences_svs"],
        rank_variants    : ["genmod_reduced_penetrance", "genmod_score_config_snvs", "genmod_score_config_svs"],
        repeat_calling   : ["str_bed"],
        repeat_annotation: ["stranger_repeat_catalog"],
    ]

    def parameterStatus = [
        workflow: [
            skip_snv_calling         : params.skip_snv_calling,
            skip_peddy               : params.skip_peddy,
            skip_phasing             : params.skip_phasing,
            skip_methylation_pileups : params.skip_methylation_pileups,
            skip_rank_variants       : params.skip_rank_variants,
            skip_repeat_calling      : params.skip_repeat_calling,
            skip_repeat_annotation   : params.skip_repeat_annotation,
            skip_snv_annotation      : params.skip_snv_annotation,
            skip_sv_calling          : params.skip_sv_calling,
            skip_sv_annotation       : params.skip_sv_annotation,
            skip_call_paralogs       : params.skip_call_paralogs,
            skip_alignment           : params.skip_alignment,
            skip_qc                  : params.skip_qc,
            skip_genome_assembly     : params.skip_genome_assembly,
        ],
        files: [
            par_regions              : params.par_regions,
            echtvar_snv_databases    : params.echtvar_snv_databases,
            svdb_sv_databases        : params.svdb_sv_databases,
            somalier_sites           : params.somalier_sites,
            vep_cache                : params.vep_cache,
            hificnv_expected_xy_cn   : params.hificnv_expected_xy_cn,
            hificnv_expected_xx_cn   : params.hificnv_expected_xx_cn,
            hificnv_excluded_regions : params.hificnv_excluded_regions,
            fasta                    : params.fasta,
            str_bed                  : params.str_bed,
            stranger_repeat_catalog  : params.stranger_repeat_catalog,
            genmod_reduced_penetrance: params.genmod_reduced_penetrance,
            genmod_score_config_snvs : params.genmod_score_config_snvs,
            genmod_score_config_svs  : params.genmod_score_config_svs,
            variant_consequences_snvs: params.variant_consequences_snvs,
            variant_consequences_svs : params.variant_consequences_svs,
            vep_plugin_files         : params.vep_plugin_files,
        ]
    ]

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters(parameterStatus, workflowSkips, workflowDependencies, fileDependencies)
    validatePacBioLicense()
    validateWorkflowCompatibility()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .ifEmpty { error "Error: No samples found in samplesheet." }
        .map { meta, reads ->
            [ meta.id, meta, reads ] // add sample as groupTuple key
        }
        .groupTuple() // group by sample
        .map {
            validateUniqueFilenamesPerSample(it)
        }
        .map { sample, metas, reads ->
            // Add number of files per sample _after_ splitting to meta
            [ sample, metas[0] + [n_files: metas.size() + metas.size() * Math.max(0, params.alignment_processes - 1), single_end:true ], reads ]
        }
        // Convert back to [ meta, reads ]
        .flatMap { _sample, meta, reads ->
            reads.collect { return [ meta, it ] }
        }
        // Add relationships to meta
        .map { meta, reads -> [ meta.family_id, meta, reads ] }
        .groupTuple()
        .map { _family, metas, reads ->
            [ addRelationshipsToMeta(metas), reads ]
        }
        .transpose()
        .set { ch_samplesheet }

        // Check that all families has at least one sample with affected phenotype if ranking is active
        validateAllFamiliesHasAffectedSamples(ch_samplesheet, params)

        // Check that there's no more than one project
        validateSingleProjectPerRun(ch_samplesheet)

        // Check that the SV calling parameters are valid
        validateSVCallingParameters()

        // Check that mothers are female, and fathers are male
        validateParentalSex(ch_samplesheet)

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//

def validateInputParameters(statusMap, workflowMap, workflowDependencies, fileDependencies) {
    validateParameterCombinations(statusMap, workflowMap, workflowDependencies, fileDependencies)

}

//
// Validate channels from input samplesheet
//
def validateUniqueFilenamesPerSample(input) {
    // Filenames needs to be unique for each sample to avoid collisions when merging
    def fileNames = input[2].collect { new File(it.toString()).name }
    if (fileNames.size() != fileNames.unique().size()) {
        error "Error: Input filenames needs to be unique for each sample."
    }
    return input
}

//
// Generate methods description for MultiQC
//
def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

def extractSoftwareFromVersions(module_yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def yamlData = yaml.load(module_yaml_file)
    // Extract all software (keys) from a module yaml
    def softwareInModule = yamlData.values().collect { it.keySet() }.flatten()
    return softwareInModule
}

def generateReferenceHTML(tool_list, description) {
    def items = tool_list
        .collect { it.trim() }
        .unique()              // samtools and bcftools share reference
        .findAll { it != "" }  // some tools does not have a reference, e.g. awk, gunzip

    if (description == 'citation') {
        return "  <p>Tools used in the workflow included: ${items.join(', ')}.</p>"
    } else if (description == 'bibliography') {
        return "  <h4>References</h4><ul><li>${items.join('</li><li>')}</li></ul>"
    }
}

def citationBibliographyText(ch_versions, references_yaml, description) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def softwareReferences = yaml.load(references_yaml.text).tool

    def unwantedReferences = ['genomic-medicine-sweden/nallo', 'Nextflow']
    // These are not collected in ch_versions but should be referenced
    def baseTools = Channel.from(['nextflow', 'nf_core', 'bioconda', 'biocontainers', 'multiqc'])

    ch_versions
        .map { module_yaml -> extractSoftwareFromVersions(module_yaml) }
        .flatten() // split multi-tool modules
        .unique()
        .filter { tool -> !unwantedReferences.contains(tool) }
        .concat(baseTools)
        .collect { tool ->
            def toolDetails = softwareReferences[tool]
            if (toolDetails == null) {
                throw new IllegalStateException("Tool: '${tool}' not found in ${references_yaml}")
            }
            return toolDetails[description]
        }
        .sort()
        .map { tools -> generateReferenceHTML(tools, description) }
}

//
// Validate  workflow skip combinations
//
def validateParameterCombinations(statusMap, workflowMap, workflowDependencies, fileDependencies) {
    // Array to store errors
    def errors = []
    // For each of the "workflow", "files"
    statusMap.each { paramsType, paramsMap ->
        paramsMap.each { param, _paramStatus ->
            if (paramsType == "files") {
                checkFileDependencies(param, fileDependencies, statusMap, workflowMap, errors)
            } else if (paramsType == "workflow") {
                checkWorkflowDependencies(param, workflowDependencies, statusMap, workflowMap, errors)
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
    def currentWorkflow = workflowMap.find { _key, mapValue -> mapValue == skip }?.key

    // If the --skip is not set, then the workflow is active, give no error
    def workflowIsActive = !statusMap["workflow"][skip]
    if(workflowIsActive) {
        return
    }

    // Get all other workflows that are required for a certain workflow
    def requiredWorkflows = combinationsMap.findAll { it.value.contains(currentWorkflow) }.keySet()
    // If no direct dependencies are found or combinationsMap does not contain the workflow, return an empty list
    if (!requiredWorkflows) {
        return []
    }
    // Collect the required --skips that are not active for the current workflow
    def dependencyString = findRequiredSkips("workflow", requiredWorkflows, statusMap, workflowMap)
        .collect { [ '--', it ].join('') }
        .join(" ")
    // If all required sets are set, give no error
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
    // Get all workflows required by a file
    def workflowThatRequiresFile = findKeysForValue(file, combinationsMap)

    workflowThatRequiresFile.each { workflow ->
        // Get the "--skip" for that workflow
        def workflowSkip = workflowMap[workflow]
        // Get the status of the "--skip", if false then workflow is active
        def WorkflowIsActive = !statusMap["workflow"][workflowSkip]
        // Get the file path
        def FilePath = statusMap["files"][file]
        // If the workflow that requires the file is active & theres no file available
        if(WorkflowIsActive && FilePath == null) {
            errors << "--$workflowSkip is NOT active, the following files are required: --$file"
        }
    }

    return errors
}

//
// Find the workflow skips that are not currently active
//
def findRequiredSkips(paramType, Set<String> requiredWorkflows, Map statusMap, Map workflowMap) {

    def requiredSkips = []

    requiredWorkflows.each { currentWorkflow ->
        // Get the skip associated with the workflow
        def skip = workflowMap[currentWorkflow]

        def workflowIsSkipped = !statusMap[paramType][skip]

        if(paramType == "workflow") {
            if(workflowIsSkipped) {
                requiredSkips << skip
            }
        }
    }
    return requiredSkips
}

def findKeysForValue(def valueToFind, Map map) {

    def keys = []

    map.each { entry ->
        def key = entry.key
        def value = entry.value

        if ((value instanceof List && value.contains(valueToFind)) || value == valueToFind) {
            keys << key
        }
    }
    return keys.isEmpty() ? null : keys
}

// Utility function to create channels from references
def createReferenceChannelFromPath(param, defaultValue = '') {
    return param ? Channel.fromPath(param, checkIfExists: true)
        .map { [ [ id: it.simpleName ], it ] }
        .collect() : defaultValue
}
// Utility function to create channels from samplesheets
def createReferenceChannelFromSamplesheet(param, schema, defaultValue = '') {
    return param ? Channel.fromList(samplesheetToList(param, schema)) : defaultValue
}

def validatePacBioLicense() {
    if (params.preset == "ONT_R10") {
        if (params.phaser.matches('hiphase')) {
            error "ERROR: The HiPhase license only permits analysis of data from PacBio."
        }
        if (params.str_caller.matches('trgt')) {
            error "ERROR: The TRGT license only permits analysis of data from PacBio."
        }
    }
}

// Genmod within RANK_VARIANTS requires affected individuals in the samplesheet.
// This is a convenience function to fail early if there are families without affected individuals.
def validateAllFamiliesHasAffectedSamples(ch_samplesheet, params) {

    if (params.skip_rank_variants) {
        return
    }

    def familiesWithPhenotypes = ch_samplesheet
        .map { meta, _reads -> [ meta.family_id, meta.phenotype ] }
        .groupTuple()

    def familiesWithoutAffected = familiesWithPhenotypes
        .filter { _family, phenotype -> !phenotype.contains(2) }

    familiesWithoutAffected
        .map { family, _phenotype -> family }
        .collect()
        .subscribe { familyList ->
            if (familyList) {
                error("ERROR: No samples in families: ${familyList.join(", ")} have affected phenotype (=2); --skip_rank_variants has to be active.")
            }
        }
}

def validateSingleProjectPerRun(ch_samplesheet) {
    ch_samplesheet
        .map { meta, _reads -> meta.project }
        .unique()
        .count()
        .map { n ->
            if ( n > 1 ) {
                error("ERROR: Only one project may be specified per run.")
            }
        }
}

def validateWorkflowCompatibility() {
    if (params.str_caller.matches('strdust') && !params.skip_repeat_annotation) {
        error "ERROR: Repeat annotation is not supported for STRdust. Run with --skip_repeat_annotation if you want to use STRdust."
    }

    if (!params.skip_sv_calling && params.sv_callers.split(',').collect { it.toLowerCase().trim() }.contains('hificnv')) {
        if (params.skip_snv_calling) {
            error "ERROR: HiFiCNV requires SNV calling to be active. Run without --skip_snv_calling if you want to use HiFiCNV."
        }
        if (!params.hificnv_expected_xy_cn || !params.hificnv_expected_xx_cn || !params.hificnv_excluded_regions) {
            error "ERROR: HiFiCNV requires expected XY and XX CN files and excluded regions to be provided. Please provide --hificnv_expected_xy_cn, --hificnv_expected_xx_cn and --hificnv_excluded_regions parameters."
        }
    }
}

def validateSVCallingParameters() {
    def sv_callers = params.sv_callers_to_merge.split(',').collect { it.toLowerCase().trim() }
    def sv_caller_priority = params.sv_callers_merge_priority.split(',').collect { it.toLowerCase().trim() }

    if (sv_callers.toSet() != sv_caller_priority.toSet()) {
        error "ERROR: The --sv_callers_merge_priority list must contain the same items as --sv_callers_to_merge (order may differ)."
    }
}

def validateParentalSex(input) {
    input
        .map { meta, _reads ->
            def sex_as_string = meta.sex == 1 ? 'male' : meta.sex == 2 ? 'female' : 'unknown'

            if ((meta.relationship == 'mother' && !isFemale(meta)) ||
                (meta.relationship == 'father' && !isMale(meta))) {
                error "ERROR: Sample ${meta.id} has been set as ${meta.relationship}, but sex is ${meta.sex} (=${sex_as_string}) in samplesheet. " +
                      "Please check the samplesheet and correct the sex or releationship."
            }
        }
}

def getParentalIds(samples, field) {
    samples.collect { it[field] }.findAll { isNonZeroNonEmpty(it) }
}

def addRelationShipsToMeta(samples) {
    // This function adds relationships to the samples based on their parental IDs.
    // We assume we are mainly interested in children, therefore if there was a grandparent, parent, and a child present,
    // the parent will be set as 'mother' or 'father' rather than as 'child'.

    def maternal_ids = getParentalIds(samples, 'maternal_id')
    def paternal_ids = getParentalIds(samples, 'paternal_id')
    def parents_ids = maternal_ids + paternal_ids
    def grandparents_ids = samples.findAll { it.id in parents_ids }.collect { it.maternal_id } +
                           samples.findAll { it.id in parents_ids }.collect { it.paternal_id }

    samples.each { sample ->
        sample.relationship = sample.id in grandparents_ids ? 'unknown' :
                              sample.id in maternal_ids ? 'mother' :
                              sample.id in paternal_ids ? 'father' :
                              isChild(sample, maternal_ids, paternal_ids) ? 'child' : 'unknown'
    }

}

def isChild(sample, maternal_ids, paternal_ids) {
    hasMother (sample, maternal_ids) ||
    hasFather (sample, paternal_ids)
}

def hasMother(sample, maternal_ids) {
    isNonZeroNonEmpty(sample.maternal_id) && (sample.maternal_id in maternal_ids)
}

def hasFather(sample, paternal_ids) {
    isNonZeroNonEmpty(sample.paternal_id) && (sample.paternal_id in paternal_ids)
}

def isFemale(sample) {
    sample.sex == 2
}

def isMale(sample) {
    sample.sex == 1
}

def boolean isNonZeroNonEmpty(value) {
    (value instanceof String && value != "" && value != "0") ||
    (value instanceof Number && value != 0)
}
