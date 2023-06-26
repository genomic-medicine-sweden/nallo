<!-- [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/fellen31/skierfe)

## Introduction

**fellen31/skierfe** is a bioinformatics analysis pipeline for long-read rare disease SV/SNV identification. Heavily influenced by best-practice pipelines such as [nf-core/nanoseq](https://github.com/nf-core/nanoseq), [nf-core/sarek](https://nf-co.re/sarek), [nf-core/raredisease](https://nf-co.re/raredisease), [PacBio Human WGS Workflow](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake), [epi2me-labs/wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) and [brentp/rare-disease-wf](https://github.com/brentp/rare-disease-wf).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

<!-- On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.
-->

## Pipeline summary

##### QC
- FastQC ([`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Aligned read QC ([`cramino`](https://github.com/wdecoster/cramino))
- Depth information ([`mosdepth`](https://github.com/brentp/mosdepth))

##### Alignment/assembly
- Align reads to reference ([`minimap2`](https://github.com/lh3/minimap2))
- Assemble (trio-binned) haploid genomes ([`hifiasm`](https://github.com/chhylp123/hifiasm))

##### Variant calling 
- Singe-sample structural variant calling and joint genotyping ([`sniffles`](https://github.com/fritzsedlazeck/Sniffles))
- Singe-sample/trio short variant calling ([`deepvariant/deeptrio`](https://github.com/google/deepvariant) or [`pepper_margin_deepvariant`](https://github.com/kishwarshafin/pepper))
    - Merge and joint genotyping of SNVs ([`GLNexus`](https://github.com/dnanexus-rnd/GLnexus))
- Assembly based variant calls ([`dipcall`](https://github.com/lh3/dipcall))
- Tandem repeat analysis if Revio/PacBio ([`TRGT`](https://github.com/PacificBiosciences/trgt/tree/main))


##### Phasing and methylation 
- Phase and haplotag reads ([`whatshap`](https://github.com/whatshap/whatshap) + [`hiphase`](https://github.com/PacificBiosciences/HiPhase))
- Methylation pileups (per haplotype), if Revio/ONT ([`modkit`](https://github.com/nanoporetech/modkit)) 

##### Annotation
- TBD
  
##### Filtering 

- TBD

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

1. Prepare a samplesheet with input data (gzipped fastq-files):

`samplesheet.csv`
```
sample,file,family_id,paternal_id,maternal_id,sex,phenotype
HG002,/path/to/HG002.fastq.gz,FAM1,HG003,HG004,1,1
HG005,/path/to/HG005.fastq.gz,FAM1,HG003,HG004,2,1
```

2. If running dipcall, download a BED file with PAR regions ([hg38](https://raw.githubusercontent.com/lh3/dipcall/master/data/hs38.PAR.bed)) and if running TRGT, download a BED file with tandem repeats ([TRGT](https://github.com/PacificBiosciences/trgt/tree/main/repeats)) matching your reference genome.

> **Note** If running dipcall, make sure chrY PAR is hard masked in reference.

3. Download the pipeline and ~~test it on a minimal dataset run with a single command~~ run:

   ```bash
   nextflow run fellen31/skierfe -r dev -profile YOURPROFILE \
     --input samplesheet.csv \
     --ped samples.ped \
     --outdir <OUTDIR> \
     --dipcall_par /path/to/hs38.PAR.bed \
     --fasta /path/to/GRCh38_no_alt_analysis_set.fasta \
     --trgt_repeats repeat_catalog_and_pathogenic.bed \
     --preset < revio/pacbio/ONT_R9/ONT_R10 >
   ```

To run in an offline environment, download the pipeline using [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use):

   ```
   nf-core download fellen31/skierfe -r dev --container singularity
   ```

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter` and `charliecloud` and which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

fellen31/skierfe was originally written by Felix Lenner.

We thank the following people for their extensive assistance in the development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, ~~please see the [contributing guidelines](.github/CONTRIBUTING.md)~~ any contribution is very welcome.

## Citations

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
