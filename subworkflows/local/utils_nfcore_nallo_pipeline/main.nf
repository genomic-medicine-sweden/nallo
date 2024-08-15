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
    aligned_read_qc  : ["mapping"],
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
// E.g., the par_regions file is required by the assembly workflow and the assembly workflow can't run without par_regions
//
def fileDependencies = [
    mapping          : ["fasta", "somalier_sites"],
    assembly         : ["fasta", "par_regions"], // The assembly workflow should be split into two - assembly and variant calling (requires ref)
    snv_calling      : ["fasta", "par_regions"],
    snv_annotation   : ["snp_db", "vep_cache", "vep_plugin_files", "reduced_penetrance", "score_config_snv", "variant_consequences_snv"],
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
        par_regions             : params.par_regions,
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
            [ sample, metas[0] + [n_files: metas.size() + metas.size() * Math.max(0, params.parallel_alignments - 1), single_end:true ], reads ]
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

        // Check that there's no more than one project
        // TODO: Try to do this in nf-schema
        ch_samplesheet
            .map { meta, reads -> meta.project }
            .unique()
            .collect()
            .filter{ it.size() == 1 }
            .ifEmpty {
                error("Only one project may be specified per run")
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

    def citation_text = [
        "MultiQC (Ewels et al. 2016)",
        "SAMtools (Danecek et al. 2021)",
    ]
    if (!params.skip_raw_read_qc) {
        citation_text = citation_text + [
            "FastQC (Andrews 2010)",
            "fcqrs",
        ]
    }
    if (!params.skip_mapping_wf) {
        if (params.parallel_alignments > 1) {
            citation_text = citation_text + [
                "splitubam",
            ]
        }
        citation_text = citation_text + [
            "SAMtools (Danecek et al. 2021)",
            "Minimap2 (Li 2018)",
            "Somalier (Pedersen et al. 2020)",
            "Sniffles2 (Smolka et al. 2024)",
        ]
        if (!params.skip_aligned_read_qc) {
            citation_text = citation_text + [
                "cramino (De Coster & Rademakers 2023)",
                "mosdepth (Pedersen & Quinlan 2018)",
            ]
        }
        if (!params.skip_call_paralogs) {
            citation_text = citation_text + [
                "paraphase",
            ]
        }
        if (!params.skip_assembly_wf) {
            if (params.hifiasm_mode == 'trio-binning') {
                citation_text = citation_text + [
                    "yak",
                ]
            }
            citation_text = citation_text + [
                "Hifiasm (Cheng et al. 2021)",
                "Gfastats (Formenti et al. 2022)",
                "dipcall (Li et al. 2018)",
                "SAMtools (Danecek et al. 2021)",
                "Minimap2 (Li 2018)",
            ]
        }
        if (!params.skip_short_variant_calling) {
            citation_text = citation_text + [
                "BEDTools (Quinlan & Hall 2010)",
                "BCFtools (Danecek et al. 2021)",
                "DeepVariant (Poplin et al. 2018)",
                "GLnexus (Yun et al. 2021)",
            ]
        }
        if (!params.skip_snv_annotation) {
            citation_text = citation_text + [
                "CADD (Rentzsch et al. 2019, Rentzsch et al. 2021)",
                "BCFtools (Danecek et al. 2021)",
                "VEP (McLaren et al. 2016)",
                "Tabix (Li 2011)",
                "Echtvar (Pedersen & de Ridder 2023)",
            ]
            if (!params.skip_rank_variants) {
                citation_text = citation_text + [
                    "Genmod (Magnusson et al. 2018)",
                    "Tabix (Li 2011)",
                ]
            }
        }
        if (!params.skip_cnv_calling) {
            citation_text = citation_text + [
                "HiFiCNV",
            ]
        }
        if (!params.skip_phasing_wf) {
            citation_text = citation_text + [
                "SAMtools (Danecek et al. 2021)",
                "cramino (De Coster & Rademakers 2023)",
            ]
            if(params.phaser == 'whatshap') {
                citation_text = citation_text + [
                    "WhatsHap (Martin et al. 2016)",
                ]
            }
            if(params.phaser == 'hiphase_sv') {
                citation_text = citation_text + [
                    "HiPhase (Holt et al. 2024)",
                ]
            }
            if(params.phaser == 'hiphase_snv') {
                citation_text = citation_text + [
                    "HiPhase (Holt et al. 2024)",
                ]
            }
            if (!params.skip_methylation_wf) {
                citation_text = citation_text + [
                    "modkit",
                    "Tabix (Li 2011)",
                ]
            }
            if (!params.skip_repeat_calling) {
                citation_text = citation_text + [
                    "TRGT (Dolzhenko et al. 2024)",
                ]
                if (!params.skip_repeat_annotation) {
                    citation_text = citation_text + [
                        "Stranger (Nilsson & Magnusson 2021)",
                    ]
                }
            }
        }
    }

    def return_text = "Tools used in the workflow included: " + citation_text.unique(false) { a, b -> a <=> b }.join(', ') - "" + "."
    return return_text
}

def toolBibliographyText() {

    reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        "<li>Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.</li>",
        "<li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.</li>",
        "<li>Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.</li>",
        "<li>Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008. doi:10.1093/gigascience/giab008</li>",
        "<li>Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842.</li>",
        "<li>Wouter De Coster, Rosa Rademakers, NanoPack2: population-scale evaluation of long-read sequencing data, Bioinformatics, Volume 39, Issue 5, May 2023, btad311, https://doi.org/10.1093/bioinformatics/btad311</li>",
        "<li>Rentzsch P, Schubach M, Shendure J, Kircher M. CADD-Splice—improving genome-wide variant effect prediction using deep learning-derived splice scores. Genome Med. 2021;13(1):31. doi:10.1186/s13073-021-00835-9</li>",
        "<li>Rentzsch P, Witten D, Cooper GM, Shendure J, Kircher M. CADD: predicting the deleteriousness of variants throughout the human genome. Nucleic Acids Research. 2019;47(D1):D886-D894. doi:10.1093/nar/gky1016</li>",
        "<li>Poplin R, Chang PC, Alexander D, et al. A universal SNP and small-indel variant caller using deep neural networks. Nat Biotechnol. 2018;36(10):983-987. doi:10.1038/nbt.4235</li>",
        "<li>Li H, Bloom JM, Farjoun Y, Fleharty M, Gauthier L, Neale B, MacArthur D (2018) A synthetic-diploid benchmark for accurate variant-calling evaluation. Nat Methods, 15:595-597. [PMID:30013044]</li>",
        "<li>Brent S Pedersen, Jeroen de Ridder, Echtvar: compressed variant representation for rapid annotation and filtering of SNPs and indels, Nucleic Acids Research, Volume 51, Issue 1, 11 January 2023, Page e3, https://doi.org/10.1093/nar/gkac931</li>",
        "<li>McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016;17(1):122. doi:10.1186/s13059-016-0974-4</li>",
        "<li>Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].</li>",
        "<li>Magnusson M, Hughes T, Glabilloy, Bitdeli Chef. genmod: Version 3.7.3. Published online November 15, 2018. doi:10.5281/ZENODO.3841142</li>",
        "<li>Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristóbal Gallardo-Alba, Alice Giani, Olivier Fedrigo, Erich D Jarvis, Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs, Bioinformatics, Volume 38, Issue 17, September 2022, Pages 4214–4216, https://doi.org/10.1093/bioinformatics/btac460</li>",
        "<li>Yun T, Li H, Chang PC, Lin MF, Carroll A, McLean CY. Accurate, scalable cohort variant calls using DeepVariant and GLnexus. Robinson P, ed. Bioinformatics. 2021;36(24):5582-5589. doi:10.1093/bioinformatics/btaa1081</li>",
        "<li>Cheng, H., Concepcion, G.T., Feng, X. et al. Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nat Methods 18, 170–175 (2021). https://doi.org/10.1038/s41592-020-01056-5</li>",
        "<li>James M Holt, Christopher T Saunders, William J Rowell, Zev Kronenberg, Aaron M Wenger, Michael Eberle, HiPhase: jointly phasing small, structural, and tandem repeat variants from HiFi sequencing, Bioinformatics, Volume 40, Issue 2, February 2024, btae042, https://doi.org/10.1093/bioinformatics/btae042</li>",
        "<li>Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191</li>",
        "<li>Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Hancock J, ed. Bioinformatics. 2018;34(5):867-868. doi:10.1093/bioinformatics/btx699</li>",
        "<li>Genome-wide profiling of highly similar paralogous genes using HiFi sequencing. Xiao Chen, Daniel Baker, Egor Dolzhenko, Joseph M Devaney, Jessica Noya, April S Berlyoung, Rhonda Brandon, Kathleen S Hruska, Lucas Lochovsky, Paul Kruszka, Scott Newman, Emily Farrow, Isabelle Thiffault, Tomi Pastinen, Dalia Kasperaviciute, Christian Gilissen, Lisenka Vissers, Alexander Hoischen, Seth Berger, Eric Vilain, Emmanuèle Délot, UCI Genomics Research to Elucidate the Genetics of Rare Diseases (UCI GREGoR) Consortium, Michael A Eberle. bioRxiv 2024.04.19.590294; doi: https://doi.org/10.1101/2024.04.19.590294</li>",
        "<li>Smolka, M., Paulin, L.F., Grochowski, C.M. et al. Detection of mosaic and population-level structural variants with Sniffles2. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-023-02024-y</li>",
        "<li>Pedersen, B.S., Bhetariya, P.J., Brown, J. et al. Somalier: rapid relatedness estimation for cancer and germline studies using efficient genome sketches. Genome Med 12, 62 (2020). https://doi.org/10.1186/s13073-020-00761-2</li>",
        "<li>Nilsson D, Magnusson M. moonso/stranger v0.7.1. Published online February 18, 2021. doi:10.5281/ZENODO.4548873</li>",
        "<li>Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics. 2011;27(5):718-719. doi:10.1093/bioinformatics/btq671</li>",
        "<li>Dolzhenko, E., English, A., Dashnow, H. et al. Characterization and visualization of tandem repeats at genome scale. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-023-02057-3</li>",
        "<li>Marcel Martin, Murray Patterson, Shilpa Garg, Sarah O Fischer, Nadia Pisanti, Gunnar W Klau, Alexander Schöenhuth, Tobias Marschall. bioRxiv 085050; doi: https://doi.org/10.1101/085050</li>",
        "<li>Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.</li>",
        "<li>Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.</li>",
        "<li>da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.</li>",
        "<li>Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.</li>",
        "<li>Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.</li>",
    ].join(' ').trim()

    return reference_text
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

