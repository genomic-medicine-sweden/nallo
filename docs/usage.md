# genomic-medicine-sweden/nallo: Usage

## Introduction

genomic-medicine-sweden/nallo is a bioinformatics analysis pipeline to analyse long-read data.

## Prerequisites

1. Install Nextflow (>=22.10.1) using the instructions [here.](https://nextflow.io/docs/latest/getstarted.html#installation)
2. Install one of the following technologies for full pipeline reproducibility: Docker, Singularity, Podman, Shifter or Charliecloud.
   > Almost all nf-core pipelines give you the option to use conda as well. However, some tools used in the nallo pipeline do not have a conda package so we do not support conda at the moment.

## Run genomic-medicine-sweden/nallo with test data

Before running the pipeline with your data, we recommend running it with the test dataset available in the `assets/test_data` folder provided with the pipeline. You do not need to download any of the data as part of it came directly with the pipeline and the other part will be fetched automatically for you when you use the test profile.

Run the following command, where YOURPROFILE is the package manager you installed on your machine. For example, `-profile test,docker` or `-profile test,singularity`:

```
nextflow run genomic-medicine-sweden/nallo \
    -profile test,<YOURPROFILE> \
    --outdir <OUTDIR>
```

> Check [nf-core/configs](https://github.com/nf-core/configs/tree/master/conf) to see if a custom config file to run nf-core pipelines already exists for your institute. If so, you can simply use `-profile test,<institute>` in your command. This enables the appropriate package manager and sets the appropriate execution settings for your machine.
> NB: The order of profiles is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

Running the command creates the following files in your working directory:

```
work                # Directory containing the Nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other Nextflow hidden files, like history of pipeline logs.
```

> [!NOTE]
> The default cpu and memory configurations used in nallo are written keeping the test profile (and dataset, which is tiny) in mind. You should override these values in configs to get it to work on larger datasets. Check the section `custom-configuration` below to know more about how to configure resources for your platform.

### Updating the pipeline

The above command downloads the pipeline from GitHub, caches it, and tests it on the test dataset. When you run the command again, it will fetch the pipeline from cache even if a more recent version of the pipeline is available. To make sure that you're running the latest version of the pipeline, update the cached version of the pipeline by including `-latest` in the command.

## Run genomic-medicine-sweden/nallo with your data

Running the pipeline involves three steps:

1. Prepare a samplesheet
2. Gather all required references
3. Supply samplesheet and references, and run the command

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

It has to be a comma-separated file with 6 columns, and a header row as shown in the examples below.
`file` can either be a gzipped-fastq file or an aligned or unalinged BAM file (BAM files will be converted to FASTQ and aligned again).
If you don't have related samples, `family_id` could be set to sample name, and `paternal_id` and `maternal_id` should be set to 0.

```console
sample,file,family_id,paternal_id,maternal_id,sex,phenotype
HG002,/path/to/HG002.fastq.gz,FAM,HG003,0,1,2
HG003,/path/to/HG003.bam,FAM,0,0,2,1
```

| Fields        | Description                                                                                                               |
| ------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `sample`      | Custom sample name, cannot contain spaces.                                                                                |
| `file`        | Absolute path to gzipped FASTQ or BAM file. File has to have the extension ".fastq.gz", .fq.gz" or ".bam".                |
| `family_id`   | "Family ID must be provided and cannot contain spaces. If no family ID is available you can use the same ID as the sample |
| `paternal_id` | Paternal ID must be provided and cannot contain spaces. If no paternal ID is available, use 0.                            |
| `maternal_id` | Maternal ID must be provided and cannot contain spaces. If no maternal ID is available, use 0.                            |
| `sex`         | Sex (0=unknown; 1=male; 2=female).                                                                                        |
| `phenotype`   | Affected status of patient (0 = missing; 1=unaffected; 2=affected).                                                       |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

The typical command for running the pipeline is as follows:

```bash
nextflow run genomic-medicine-sweden/nallo -profile docker \
    --input samplesheet.csv \
    --preset <revio/pacbio/ONT_R10> \
    --outdir <OUTDIR> \
    --fasta <reference.fasta> \
    --skip_assembly_wf \
    --skip_repeat_wf \
    --skip_snv_annotation \
    --skip_cnv_calling
```

## Presets

This pipeline comes with three different presets that can be set with the `--preset` parameter, "revio", "pacbio" and "ONT_R10" (defaults to "revio").
By default, `--skip_assembly_wf` and `--skip_repeat_wf` is set to true for `ONT_R10` and `--skip_methylation_wf` is set to true for "pacbio", which means these workflows are not run.

## Reference files and parameters

The typical command example above requires no additional files except the reference genome.
Nallo has the ability to skip certain parts of the pipeline, for example `--skip_repeat_wf`.
Some workflows require additional files:

- If running without `--skip_assembly_wf`, download a BED file with PAR regions ([hg38](https://raw.githubusercontent.com/lh3/dipcall/master/data/hs38.PAR.bed)) to supply with `--dipcall_par`.

> [!NOTE]
> Make sure chrY PAR is hard masked in reference.

- If running without `--skip_repeat_calling`, download a BED file with tandem repeats ([TRGT](https://github.com/PacificBiosciences/trgt/tree/main/repeats)) matching your reference genome to supply with `--trgt_repeats`.

- If running without `--skip_repeat_annotation`, download a json variant catalog, (e.g. [variant_catalog_grch38.json](https://github.com/Clinical-Genomics/stranger/raw/main/stranger/resources/variant_catalog_grch38.json)) matching your reference genome to supply with `--variant_catalog`.

- If running without `--skip_snv_annotation`, download [VEP cache](https://ftp.ensembl.org/pub/release-110/variation/vep/homo_sapiens_vep_110_GRCh38.tar.gz) to supply with `--vep_cache` and prepare a samplesheet with annotation databases ([`echtvar encode`](https://github.com/brentp/echtvar)) to supply with `--snp_db`:
  - If your samplesheet contains at least one affected sample (phenotype = 2), `--reduced_penetrance` (Used by GENMOD while modeling the variants. Contains a list of loci that show [reduced penetrance](https://medlineplus.gov/genetics/understanding/inheritance/penetranceexpressivity/) in people. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/reduced_penetrance.tsv)), `--score_config_snv` (Used by GENMOD for ranking the variants. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/rank_model_snv.ini)) and `--variant_consequences_snv` (File containing list of SO terms listed in the order of severity from most severe to lease severe for annotating genomic and mitochondrial SNVs. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/variant_consequences_v2.txt). You can learn more about these terms [here](https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html)) is also required.
    `snp_dbs.csv`

```
sample,file
gnomad,/path/to/gnomad.v3.1.2.echtvar.popmax.v2.zip
cadd,/path/to/cadd.v1.6.hg38.zip
```

- If running without `--skip_cnv_calling`, expected CN regions for your reference genome can be downloaded from [HiFiCNV GitHub](https://github.com/PacificBiosciences/HiFiCNV/tree/main/data) to supply with `--hificnv_xy`, `--hificnv_xx` (expected_cn) and `--hificnv_exclude` (excluded_regions).

- If you want to include extra samples for mili-sample calling of SVs - prepare a samplesheet with .snf files from Sniffles to supply with `--extra_snfs`:

`extra_snfs.csv`

```
sample,file
HG01123,/path/to/HG01123_sniffles.snf
HG01124,/path/to/HG01124_sniffles.snf
```

- If running without `--skip_call_paralogs`, the reference genome needs to be hg38

- If running without `--skip_mapping_wf`, a VCF of known polymorphic sites (e.g. [sites.hg38.vcg.gz](https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz)) needs to be supplied with `--somalier_sites`, from which sex will be inferred if possible.

#### Highlighted parameters:

- You can choose to limit SNV calling to regions in BED file (`--bed`).

- By default SNV-calling is split into 13 parallel processes, limit this by setting `--parallel_snv` to a different number.

- By default the pipeline does not perform parallel alignment, but this can be set by setting `--split_fastq` to split the input and alignment into N files/processes.

All parameters are listed below:

## Workflow skip options

Options to skip various steps within the workflow

| Parameter                    | Description                                | Type      | Default | Required | Hidden |
| ---------------------------- | ------------------------------------------ | --------- | ------- | -------- | ------ |
| `skip_aligned_read_qc`       | Skip aligned read QC                       | `boolean` | `False` |          |        |
| `skip_raw_read_qc`           | Skip raw read QC                           | `boolean` | `False` |          |        |
| `skip_short_variant_calling` | Skip short variant calling                 | `boolean` | `False` |          |        |
| `skip_assembly_wf`           | Skip assembly and downstream processes     | `boolean` | `False` |          |        |
| `skip_mapping_wf`            | Skip read mapping and downstream processes | `boolean` | `False` |          |        |
| `skip_methylation_wf`        | Skip methylation workflow                  | `boolean` | `False` |          |        |
| `skip_repeat_calling         | Skip repeat calling workflow               | `boolean` | `False` |          |        |
| `skip_repeat_annotation`     | Skip repeat annotation workflow            | `boolean` | `False` |          |        |
| `skip_phasing_wf`            | Skip phasing workflow                      | `boolean` | `False` |          |        |
| `skip_snv_annotation`        | Skip SNV annotation                        | `boolean` | `False` |          |        |
| `skip_cnv_calling`           | Skip CNV workflow                          | `boolean` | `False` |          |        |
| `skip_call_paralogs`         | Skip call paralogs (Paraphase)             | `boolean` | `False` |          |        |

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter       | Description                                                                                                                                                                                                                                                                                                                                                                                  | Type     | Default | Required | Hidden |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- | ------- | -------- | ------ |
| `input`         | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details> | `string` |         | True     |        |
| `outdir`        | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.                                                                                                                                                                                                                                                                     | `string` |         | True     |        |
| `email`         | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>                                  | `string` |         |          |        |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified.                                                                                                                                                                                                                                                                                                  | `string` |         |          |        |

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter         | Description                                                                                                                                                                                                                                                                                                                                                                                                                   | Type      | Default | Required | Hidden |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `fasta`           | Reference genome                                                                                                                                                                                                                                                                                                                                                                                                              | `string`  |         |          |        |
| `genome`          | Name of iGenomes reference. <details><summary>Help</summary><small>If using a reference genome configured in the pipeline using iGenomes, use this parameter to give the ID for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`. <br><br>See the [nf-core website docs](https://nf-co.re/usage/reference_genomes) for more details.</small></details> | `string`  |         |          |        |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>                                                                                                                                             | `boolean` | True    |          | True   |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter                    | Description                                                                                                                                                                                                                                                                                                                                                                                       | Type     | Default                                                  | Required | Hidden |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- | -------------------------------------------------------- | -------- | ------ |
| `custom_config_version`      | Git commit id for Institutional configs.                                                                                                                                                                                                                                                                                                                                                          | `string` | master                                                   |          | True   |
| `custom_config_base`         | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details> | `string` | https://raw.githubusercontent.com/nf-core/configs/master |          | True   |
| `config_profile_name`        | Institutional config name.                                                                                                                                                                                                                                                                                                                                                                        | `string` |                                                          |          | True   |
| `config_profile_description` | Institutional config description.                                                                                                                                                                                                                                                                                                                                                                 | `string` |                                                          |          | True   |
| `config_profile_contact`     | Institutional config contact information.                                                                                                                                                                                                                                                                                                                                                         | `string` |                                                          |          | True   |
| `config_profile_url`         | Institutional config URL link.                                                                                                                                                                                                                                                                                                                                                                    | `string` |                                                          |          | True   |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter    | Description                                                                                                                                                                                                                                                                 | Type      | Default | Required | Hidden |
| ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `max_cpus`   | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>                                      | `integer` | 16      |          | True   |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details> | `string`  | 128.GB  |          | True   |
| `max_time`   | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>        | `string`  | 240.h   |          | True   |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter                      | Description                                                                                                                                                                                                                                                                                                                                                                                                  | Type      | Default | Required | Hidden |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------- | ------- | -------- | ------ |
| `help`                         | Display help text.                                                                                                                                                                                                                                                                                                                                                                                           | `boolean` |         |          | True   |
| `version`                      | Display version and exit.                                                                                                                                                                                                                                                                                                                                                                                    | `boolean` |         |          | True   |
| `publish_dir_mode`             | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details> | `string`  | copy    |          | True   |
| `email_on_fail`                | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>                                                                                                                                                  | `string`  |         |          | True   |
| `plaintext_email`              | Send plain-text email instead of HTML.                                                                                                                                                                                                                                                                                                                                                                       | `boolean` |         |          | True   |
| `max_multiqc_email_size`       | File size limit when attaching MultiQC reports to summary emails.                                                                                                                                                                                                                                                                                                                                            | `string`  | 25.MB   |          | True   |
| `monochrome_logs`              | Do not use coloured log outputs.                                                                                                                                                                                                                                                                                                                                                                             | `boolean` |         |          | True   |
| `hook_url`                     | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>                                                                                                                                                                                                                       | `string`  |         |          | True   |
| `multiqc_config`               | Custom config file to supply to MultiQC.                                                                                                                                                                                                                                                                                                                                                                     | `string`  |         |          | True   |
| `multiqc_logo`                 | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file                                                                                                                                                                                                                                                                                                                 | `string`  |         |          | True   |
| `multiqc_methods_description`  | Custom MultiQC yaml file containing HTML including a methods description.                                                                                                                                                                                                                                                                                                                                    | `string`  |         |          |        |
| `validate_params`              | Boolean whether to validate parameters against the schema at runtime                                                                                                                                                                                                                                                                                                                                         | `boolean` | True    |          | True   |
| `validationShowHiddenParams`   | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>                                                                                                                    | `boolean` |         |          | True   |
| `validationSkipDuplicateCheck` | nf-validation related parameter                                                                                                                                                                                                                                                                                                                                                                              | `boolean` |         |          | True   |
| `validationS3PathCheck`        | Boolean whether to validate validate AWS S3 paths                                                                                                                                                                                                                                                                                                                                                            | `boolean` |         |          | True   |
| `monochromeLogs`               | Boolean whether to color nf-validation logs                                                                                                                                                                                                                                                                                                                                                                  | `boolean` |         |          | True   |

## Workflow options

| Parameter        | Description                                | Type      | Default     | Required | Hidden |
| ---------------- | ------------------------------------------ | --------- | ----------- | -------- | ------ |
| `preset`         | Choose a preset depending on data type     | `string`  | revio       | True     |        |
| `variant_caller` | Choose variant caller                      | `string`  | deepvariant |          |        |
| `phaser`         | Choose phasing software                    | `string`  | whatshap    |          |        |
| `hifiasm_mode`   | Run hifiasm in hifi-only or hifi-trio mode | `string`  | hifi-only   |          |        |
| `split_fastq`    | Split alignment into n jobs                | `integer` | 0           |          |        |
| `parallel_snv`   | Split SNV calling into n chunks            | `integer` | 13          |          |        |

## Extra file inputs

Different processes may need extra input files

| Parameter                          | Description                                                                                                                                                                                                                                                               | Type      | Default | Required | Hidden |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `dipcall_par`                      | Provide a bed file of chrX PAR regions for dipcall                                                                                                                                                                                                                        | `string`  |         |          |        |
| `extra_snfs`                       | Extra input files for Sniffles                                                                                                                                                                                                                                            | `string`  |         |          |        |
| `tandem_repeats`                   | Tandem repeat BED-file for sniffles                                                                                                                                                                                                                                       | `string`  |         |          |        |
| `trgt_repeats`                     | BED-file for repeats to be genotyped                                                                                                                                                                                                                                      | `string`  |         |          |        |
| `snp_db`                           | Extra echtvar-databases to annotate SNVs with                                                                                                                                                                                                                             | `string`  |         |          |        |
| `vep_cache`                        | Path to directory of vep_cache                                                                                                                                                                                                                                            | `string`  |         |          |        |
| `bed`                              | BED file with regions of interest                                                                                                                                                                                                                                         | `string`  |         |          |        |
| `hificnv_xy`                       |                                                                                                                                                                                                                                                                           | `string`  |         |          |        |
| `hificnv_xx`                       |                                                                                                                                                                                                                                                                           | `string`  |         |          |        |
| `hificnv_exclude`                  | HiFiCNV BED file specifying regions to exclude                                                                                                                                                                                                                            | `string`  |         |          |        |
| `somalier_sites`                   | A VCF of known polymorphic sites                                                                                                                                                                                                                                          | `string`  |         |          |        |
| `variant_catalog`                  | Variant catalog json-file for Stranger                                                                                                                                                                                                                                    | `string`  |         |          |        |
| `validationFailUnrecognisedParams` | Validation of parameters fails when an unrecognised parameter is found. <details><summary>Help</summary><small>By default, when an unrecognised parameter is found, it returns a warning.</small></details>                                                               | `boolean` |         |          | True   |
| `validationLenientMode`            | Validation of parameters in lenient more. <details><summary>Help</summary><small>Allows string values that are parseable as numbers or booleans. For further information see [JSONSchema docs](https://github.com/everit-org/json-schema#lenient-mode).</small></details> | `boolean` |         |          | True   |

### Updating the pipeline

```bash
nextflow pull genomic-medicine-sweden/nallo
```

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [genomic-medicine-sweden/nallo releases page](https://github.com/genomic-medicine-sweden/nallo/releases) and find the latest pipeline version - numeric only (eg. `0.1.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 0.1.0`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline without internet access

The pipeline and container images can be downloaded using [nf-core tools](https://nf-co.re/docs/usage/offline). For running offline, you of course have to make all the reference data available locally, and specify `--fasta`, etc., see [above](#reference-files-and-parameters).

Contrary to the paragraph about [Nextflow](https://nf-co.re/docs/usage/offline#nextflow) on the page linked above, it is not possible to use the "-all" packaged version of Nextflow for this pipeline. The online version of Nextflow is necessary to support the necessary nextflow plugins. Download instead the file called just `nextflow`. Nextflow will download its dependencies when it is run. Additionally, you need to download the nf-validation plugin explicitly:

```
./nextflow plugin install nf-validation
```

Now you can transfer the `nextflow` binary as well as its directory `$HOME/.nextflow` to the system without Internet access, and use it there. It is necessary to use an explicit version of `nf-validation` offline, or Nextflow will check for the most recent version online. Find the version of nf-validation you downloaded in `$HOME/.nextflow/plugins`, then specify this version for `nf-validation` in your configuration file:

```
plugins {
        // Set the plugin version explicitly, otherwise nextflow will look for the newest version online.
        id 'nf-validation@1.1.3'
}
```

This should go in your Nextflow confgiguration file, specified with `-c <YOURCONFIG>` when running the pipeline.
