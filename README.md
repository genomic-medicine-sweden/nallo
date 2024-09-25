[![GitHub Actions CI Status](https://github.com/genomic-medicine-sweden/nallo/actions/workflows/ci.yml/badge.svg)](https://github.com/genomic-medicine-sweden/nallo/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/genomic-medicine-sweden/nallo/actions/workflows/linting.yml/badge.svg)](https://github.com/genomic-medicine-sweden/nallo/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13748210.svg)](https://doi.org/10.5281/zenodo.13748210)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/genomic-medicine-sweden/nallo)

## Introduction

**genomic-medicine-sweden/nallo** is a bioinformatics analysis pipeline for long-reads from both PacBio and (targeted) ONT-data, focused on rare-disease. Heavily influenced by best-practice pipelines such as [nf-core/sarek](https://nf-co.re/sarek), [nf-core/raredisease](https://nf-co.re/raredisease), [nf-core/nanoseq](https://github.com/nf-core/nanoseq), [PacBio Human WGS Workflow](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake), [epi2me-labs/wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) and [brentp/rare-disease-wf](https://github.com/brentp/rare-disease-wf).

## Overview

<picture align="center">
    <img alt="genomic-medicine-sweden/nallo workflow" src="docs/images/nallo_metromap.png">
  </picture>

## Pipeline summary

##### QC

- Read QC with [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [cramino](https://github.com/wdecoster/cramino) and [mosdepth](https://github.com/brentp/mosdepth)

##### Alignment & assembly

- Align reads to reference with [minimap2](https://github.com/lh3/minimap2)
- Assemble (trio-binned) haploid genomes with [hifiasm](https://github.com/chhylp123/hifiasm) (HiFi only)

##### Variant calling

- Call SNVs & joint genotyping with [deepvariant](https://github.com/google/deepvariant) and [GLNexus](https://github.com/dnanexus-rnd/GLnexus)
- Call SVs with [Severus](https://github.com/KolmogorovLab/Severus) or [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
- Call CNVs with [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV)
- Call tandem repeats with [TRGT](https://github.com/PacificBiosciences/trgt/tree/main) (HiFi only)
- Call paralogous genes with [Paraphase](https://github.com/PacificBiosciences/paraphase)
- Call variants from assembly with [dipcall](https://github.com/lh3/dipcall) (HiFi only)

##### Phasing and methylation

- Phase and haplotag reads with [LongPhase](https://github.com/twolinin/longphase), [whatshap](https://github.com/whatshap/whatshap) or [HiPhase](https://github.com/PacificBiosciences/HiPhase)
- Create methylation pileups with [modkit](https://github.com/nanoporetech/modkit)

##### Annotation

- Annotate SNVs and INDELs with databases of choice, i.e. [gnomAD](https://gnomad.broadinstitute.org), [CADD](https://cadd.gs.washington.edu) etc. with [echtvar](https://github.com/brentp/echtvar) and [VEP](https://github.com/Ensembl/ensembl-vep)
- Annotate repeat expansions with [stranger](https://github.com/Clinical-Genomics/stranger)

##### Ranking

- Rank SNVs with [GENMOD](https://github.com/Clinical-Genomics/genmod)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet with input data:

`samplesheet.csv`

```
project,sample,file,family_id,paternal_id,maternal_id,sex,phenotype
NIST,HG002,/path/to/HG002.fastq.gz,FAM1,HG003,HG004,1,2
NIST,HG005,/path/to/HG005.bam,FAM1,HG003,HG004,2,1
```

Supply a reference genome with `--fasta` and choose a matching `--preset` for your data (`revio`, `pacbio`, `ONT_R10`). Now, you can run the pipeline using:

```bash
nextflow run genomic-medicine-sweden/nallo \
    -profile <docker/singularity/.../institute> \
    --input samplesheet.csv \
    --preset <revio/pacbio/ONT_R10> \
    --fasta <reference.fasta> \
    --outdir <OUTDIR>
```

For more details and further functionality, please refer to the [usage documentation](https://github.com/genomic-medicine-sweden/nallo/blob/dev/docs/usage.md).

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

To run in an offline environment, download the pipeline and singularity images using [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use):

```
nf-core download genomic-medicine-sweden/nallo
```

## Credits

genomic-medicine-sweden/nallo was originally written by Felix Lenner.

We thank the following people for their extensive assistance in the development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use genomic-medicine-sweden/nallo for your analysis, please cite it using the following doi: [10.5281/zenodo.13748210](https://doi.org/10.5281/zenodo.13748210).

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
