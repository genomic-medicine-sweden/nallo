# genomic-medicine-sweden/nallo: Usage

## Introduction

genomic-medicine-sweden/nallo is a bioinformatics analysis pipeline to analyse long-read data.

## Prerequisites

1. Install Nextflow (>=24.04.2) using the instructions [here.](https://nextflow.io/docs/latest/getstarted.html#installation)
2. Install one of the following technologies for full pipeline reproducibility: Docker, Singularity, Podman, Shifter or Charliecloud.
   > Almost all nf-core pipelines give you the option to use conda as well. However, some tools used in genomic-medicine-sweden/nallo do not have a conda package so we do not support conda at the moment.

## Getting started

Before running the pipeline with your data, we recommend running it with the test profile. You do not need to download any of the data as it will be fetched automatically for you when you use the test profile.

Run the following command, where YOURPROFILE is the package manager you installed on your machine. For example, `-profile test,docker` or `-profile test,singularity`

```
nextflow run genomic-medicine-sweden/nallo \
    -profile test,<YOURPROFILE> \
    --outdir <OUTDIR>
```

> Check [nf-core/configs](https://github.com/nf-core/configs/tree/master/conf) to see if a custom config file to run nf-core pipelines already exists for your institute. If so, you can simply use `-profile test,<institute>` in your command. This enables the appropriate package manager and sets the appropriate execution settings for your machine.
> NB: The order of profiles is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

Running the command creates the following files in your working directory

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

## Running genomic-medicine-sweden/nallo with your data

Running the pipeline on real data involves three steps:

1. Preparing a samplesheet with your data
2. Gather required files and references
3. Supply samplesheet, refeferences and files and run the pipeline

## Samplesheet

First, you will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

It has to be a comma-separated file with 7 columns, and a header row as shown in the example below:

```console
project,sample,file,family_id,paternal_id,maternal_id,sex,phenotype
testrun,HG002,/path/to/HG002.fastq.gz,FAM,HG003,0,1,2
testrun,HG003,/path/to/HG003.bam,FAM,0,0,2,1
```

| Fields        | Description                                                                                                                       |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `project`     | Project name must be provided and cannot contain spaces, needs to be the same for all samples."                                   |
| `sample`      | Custom sample name, cannot contain spaces.                                                                                        |
| `file`        | Absolute path to gzipped FASTQ or BAM file. File has to have the extension ".fastq.gz", .fq.gz" or ".bam".                        |
| `family_id`   | Family ID must be provided and cannot contain spaces. If no family ID is available use the same ID as sample.                     |
| `paternal_id` | Paternal ID must be provided and cannot contain spaces. If no paternal ID is available, use 0.                                    |
| `maternal_id` | Maternal ID must be provided and cannot contain spaces. If no maternal ID is available, use 0.                                    |
| `sex`         | Sex must be provided as 0, 1 or 2 (0=unknown; 1=male; 2=female). If sex is unknown it will be assigned automatically if possible. |
| `phenotype`   | Affected status of patient (0 = missing; 1=unaffected; 2=affected).                                                               |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Preset

This pipeline comes with three different presets that should be set with the `--preset` parameter

- `revio` (default)
- `pacbio`
- `ONT_R10`

`--skip_assembly_wf` and `--skip_repeat_wf` will be set to true for `ONT_R10` and `--skip_methylation_wf` will be set to true for `pacbio`, meaning these subworkflows are not run.

## Subworkflows

As indicated above, this pipeline is divided into multiple subworkflows, each with its own input requirements and outputs. By default, all subworklows are active, and thus all mandatory input files are required.

The only parameter mandatory for all subworkflows is the `--input` and `--outdir` parameters, all other parameters are determined by the active subworkflows. If you would run `nextflow run genomic-medicine-sweden/nallo -profile docker --outdir results --input samplesheet.csv`

```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  --skip_assembly_wf is NOT active, the following files are required: --dipcall_par
  --skip_snv_annotation is NOT active, the following files are required: --snp_db
  --skip_mapping_wf is NOT active, the following files are required: --somalier_sites
  --skip_snv_annotation is NOT active, the following files are required: --vep_cache
  ...
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

The pipeline will try to guide you through which files are required, but a thorough description is provided below.

Additionally, if you want to skip a subworkflow, you will need to explicitly state to skip all subworklow that relies on it. For example, `nextflow run genomic-medicine-sweden/nallo -profile docker --outdir results --input samplesheet.csv --skip_mapping_wf` will tell you

```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  --skip_mapping_wf is active, the pipeline has to be run with: --skip_qc --skip_assembly_wf --skip_call_paralogs --skip_short_variant_calling --skip_snv_annotation --skip_cnv_calling --skip_phasing_wf --skip_rank_variants --skip_repeat_calling --skip_repeat_annotation --skip_methylation_wf
  ...
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

Because almost all other subworkflows relies on the mapping subworkflow.

## Reference files and parameters

As descibed above, the files required depend on the active subworkflows. All parameters are listed [here](parameters.md), but the most useful parameters needed to run the pipeline described in more detail below.

### Mapping (`--skip_mapping_wf`)

The majority of subworkflows depend on the mapping (alignment) subworkflow which requires `--fasta` and `--somalier_sites`.

| Parameter        | Description                                                                                                                                                                              |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `fasta`          | Reference genome, either gzipped or uncompressed FASTA (e.g. [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)) |
| `somalier_sites` | A VCF of known polymorphic sites (e.g. [sites.hg38.vcg.gz](https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz)), from which sex will be inferred if possible.            |

### QC (`--skip_qc`)

This subworkflow depends on the mapping subworkflow, but requires no additional files.

### Assembly (`--skip_assembly_wf`)

This subworkflow contains both genome assembly and assembly variant calling. The assemblyt variant calling needs the sex of samples and for samples with unknown sex this is inferred from aligned reads, therefore it depends on the mapping subworkflow.

It requires a BED file with PAR regions.

| Parameter     | Description                                                                                                                        |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `par_regions` | A BED file with PAR regions (e.g. [GRCh38_PAR.bed](https://storage.googleapis.com/deepvariant/case-study-testdata/GRCh38_PAR.bed)) |

> [!NOTE]
> Make sure chrY PAR is hard masked in reference genome you are using.

### Call paralogs (`--skip_call_paralogs`)

This subworkflow depends on the mapping subworkflow, but requires no additional files.

> [!NOTE]
> Only GRCh38 is supported.

### Short variant calling (`--skip_short_variant_calling`)

This subworkflow depends on the mapping subworkflow, and required the same PAR regions file as the assembly workflow.

| Parameter     | Description                                                                                                                       |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `par_regions` | A BED file with PAR regions (e.g. [GRCh38_PAR.bed](ttps://storage.googleapis.com/deepvariant/case-study-testdata/GRCh38_PAR.bed)) |

### CNV calling (`--skip_cnv_calling`)

This subworkflow depends on the mapping and short variant calling subworkflows, and requires the following additional files:

| Parameter         | Description                                                                                                                                                                                     |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `hificnv_xy`      | expected XY copy number regions for your reference genome (e.g. [expected_cn.hg38.XY.bed](https://github.com/PacificBiosciences/HiFiCNV/raw/main/data/expected_cn/expected_cn.hg38.XY.bed))     |
| `hificnv_xx`      | expected XX copy number regions for your reference genome (e.g. [expected_cn.hg38.XX.bed](https://github.com/PacificBiosciences/HiFiCNV/raw/main/data/expected_cn/expected_cn.hg38.XX.bed))     |
| `hificnv_exclude` | BED file specifying regions to exclude (e.g. [cnv.excluded_regions.hg38.bed.gz](https://github.com/PacificBiosciences/HiFiCNV/raw/main/data/excluded_regions/cnv.excluded_regions.hg38.bed.gz)) |

### Phasing (`--skip_phasing_wf`)

This subworkflow phases variants and haplotags aligned BAM files, and such relies on the mapping and short variant calling subworkflows, but requires no additional files.

### Methylation (`--skip_methylation_wf`)

This subworkflow relies on mapping, short variant calling and phasing subworkflows, but requires no additional files.

### Repeat calling (`--skip_repeat_calling`)

This subworkflow requires haplotagged BAM files, and such relies on the mapping, short variant calling and phasing subworkflows, and requires the following additional files:

| Parameter      | Description                                                                                                                                                                                   |
| -------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `trgt_repeats` | a BED file with tandem repeats matching your reference genome (e.g. [pathogenic_repeats.hg38.bed](https://github.com/PacificBiosciences/trgt/raw/main/repeats/pathogenic_repeats.hg38.bed)>)) |

### Repeat annotation (`--skip_repeat_annotation`)

This subworkflow relies on the mapping, short variant calling, phasing and repeat calling subworkflows, and requires the following additional files:

| Parameter         | Description                                                                                                                                                                           |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `variant_catalog` | a variant catalog matching your reference (e.g. [variant_catalog_grch38.json](https://github.com/Clinical-Genomics/stranger/raw/main/stranger/resources/variant_catalog_grch38.json)) |

### SNV annotation (`--skip_snv_annotation`)

This subworkflow relies on the mapping and short variant calling, and requires the following additional files:

<!-- TODO: score_config_snv, reduced_penetrance and variant_consequences_snv should link to real examples -->

| Parameter                  | Description                                                                                                                                                                                                                                                                                                                                                                   |
| -------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `vep_cache`                | VEP cache matching your reference genome, either as a `.tar.gz` archive or path to a directory (e.g. [homo_sapiens_vep_110_GRCh38.tar.gz](https://ftp.ensembl.org/pub/release-110/variation/vep/homo_sapiens_vep_110_GRCh38.tar.gz))                                                                                                                                          |
| `vep_plugins` <sup>1</sup> | A csv file with VEP plugin files, pLI and LoFtool are required. Example provided below.                                                                                                                                                                                                                                                                                       |
| `snp_db` <sup>2</sup>      |  A csv file with annotation databases from ([`echtvar encode`](https://github.com/brentp/echtvar))                                                                                                                                                                                                                                                                            |
| `variant_consequences_snv` | A list of SO terms listed in the order of severity from most severe to lease severe for annotating genomic and mitochondrial SNVs. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/variant_consequences_v2.txt). You can learn more about these terms [here](https://ensembl.org/info/genome/variation/prediction/predicted_data.html) |

<sup>1</sup> Example file for input with `--vep_plugins`

```
vep_files
https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/nallo/reference/vep_plugins/spliceai_21_scores_raw_indel_-v1.3-.vcf.gz.tbi
https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/nallo/reference/vep_plugins/spliceai_21_scores_raw_indel_-v1.3-.vcf.gz
https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/nallo/reference/vep_plugins/spliceai_21_scores_raw_snv_-v1.3-.vcf.gz.tbi
https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/nallo/reference/vep_plugins/pLI_values.txt
https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/nallo/reference/vep_plugins/LoFtool_scores.txt
```

<sup>2</sup> Example file for input with `--snp_db`:

```
sample,file
gnomad,/path/to/gnomad.v3.1.2.echtvar.popmax.v2.zip
cadd,/path/to/cadd.v1.6.hg38.zip
```

> [!WARNING]
> Generating an echtvar database from a VCF-file is a fairly straightforward process described on the [echtvar GitHub](https://github.com/brentp/echtvar). However, the pre-made `gnomad.v3.1.2.echtvar.v2.zip` provided by them results in malformed INFO lines that are not compatible with genmod (run in the subsequent ranking subworkflow).
>
> For a very small test database that only overlaps the coordinates of the pipeline test data set, you could use [`cadd.v1.6.hg38.test_data.zip`](https://github.com/genomic-medicine-sweden/test-datasets/raw/refs/heads/nallo/reference/cadd.v1.6.hg38.test_data.zip) to get started.

> [!NOTE]
> Optionally, to calcuate CADD scores for small indels, supply a path to a folder containing cadd annotations with `--cadd_resources` and prescored indels with `--cadd_prescored`. Equivalent of the `data/annotations/` and `data/prescored/` folders described [here](https://github.com/kircherlab/CADD-scripts/#manual-installation). CADD scores for SNVs can be annotated through echvtvar and `--snp_db`.

### SV annotation (`--skip_sv_annotation`)

This subworkflow relies on the mapping subworkflow, and requires the following additional files:

| Parameter               | Description                                                                                                                                                                                                                                            |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `svdb_dbs` <sup>1</sup> | Csv file with databases used for structural variant annotation in vcf format. <details><summary>Help</summary><small>Path to comma-separated file containing information about the databases used for structural variant annotation.</small></details> |

<sup>1</sup> Example file for input with `--svdb_dbs`:

```
filename,in_freq_info_key,in_allele_count_info_key,out_freq_info_key,out_allele_count_info_key
https://github.com/genomic-medicine-sweden/test-datasets/raw/b9ff54b59cdd39df5b6e278a30b08d94075a644c/reference/colorsdb.test_data.vcf.gz,AF,AC,colorsdb_af,colorsdb_ac
```

These databases could for example come from [CoLoRSdb](https://zenodo.org/records/13145123).

### Rank variants (`--skip_rank_variants`)

This subworkflow ranks SNVs, and relies on the mapping, short variant calling and SNV annotation subworkflows, and requires the following additional files:

| Parameter            | Description                                                                                                                                                                                                                                                 |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `score_config_snv`   |  Used by GENMOD when ranking variants. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/rank_model_snv.ini).                                                                                                          |
| `reduced_penetrance` | A list of loci that show [reduced penetrance](https://medlineplus.gov/genetics/understanding/inheritance/penetranceexpressivity/) in people. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/reduced_penetrance.tsv) |

### Other highlighted parameters

- Limit SNV calling to regions in BED file (`--bed`).
- By default SNV-calling is split into 13 parallel processes, this speeds up the variant calling significantly. Limit this by setting `--parallel_snv` to a different number.
- By default the pipeline splits the input files into eight pieces, performs parallel alignment and then merges the files. This can be changed to a different number with `--parallel_alignments`, or turned off by supplying a value of 1. Parallel alignment comes with some additional overhead, but can speed up the pipeline significantly.

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [genomic-medicine-sweden/nallo releases page](https://github.com/genomic-medicine-sweden/nallo/releases) and find the latest pipeline version - numeric only (eg. `0.2.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 0.2.0`. Of course, you can switch to another version by changing the number after the `-r` flag.

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
