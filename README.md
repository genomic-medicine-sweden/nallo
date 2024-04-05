[![GitHub Actions CI Status](https://github.com/genomic-medicine-sweden/skierfe/actions/workflows/ci.yml/badge.svg)](https://github.com/genomic-medicine-sweden/skierfe/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/genomic-medicine-sweden/skierfe/actions/workflows/linting.yml/badge.svg)](https://github.com/genomic-medicine-sweden/skierfe/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/genomic-medicine-sweden/skierfe)

## Introduction

**genomic-medicine-sweden/skierfe** is a bioinformatics analysis pipeline for long-read rare disease SV/SNV identification using both PacBio and (targeted) ONT-data. Heavily influenced by best-practice pipelines such as [nf-core/nanoseq](https://github.com/nf-core/nanoseq), [nf-core/sarek](https://nf-co.re/sarek), [nf-core/raredisease](https://nf-co.re/raredisease), [PacBio Human WGS Workflow](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake), [epi2me-labs/wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) and [brentp/rare-disease-wf](https://github.com/brentp/rare-disease-wf).

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

<!-- On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.
-->

## Pipeline summary

##### QC

- FastQC ([`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Aligned read QC ([`cramino`](https://github.com/wdecoster/cramino))
- Depth information ([`mosdepth`](https://github.com/brentp/mosdepth))

##### Alignment & assembly

- Align reads to reference ([`minimap2`](https://github.com/lh3/minimap2))
- Assemble (trio-binned) haploid genomes (HiFi only) ([`hifiasm`](https://github.com/chhylp123/hifiasm))

##### Variant calling

- Short variant calling & joint genotyping of SNVs ([`deepvariant`](https://github.com/google/deepvariant) + [`GLNexus`](https://github.com/dnanexus-rnd/GLnexus))
- SV calling and joint genotyping ([`sniffles2`](https://github.com/fritzsedlazeck/Sniffles))
- Tandem repeats ([`TRGT`](https://github.com/PacificBiosciences/trgt/tree/main))
- Assembly based variant calls (HiFi only) ([`dipcall`](https://github.com/lh3/dipcall))
- CNV-calling (HiFi only) ([`HiFiCNV`](https://github.com/PacificBiosciences/HiFiCNV))

##### Phasing and methylation

- Phase and haplotag reads ([`whatshap`](https://github.com/whatshap/whatshap) + [`hiphase`](https://github.com/PacificBiosciences/HiPhase))
- Methylation pileups (Revio/ONT) ([`modkit`](https://github.com/nanoporetech/modkit))

##### Annotation - SNV

1. Annotate variants with database(s) of choice, i.e. [gnomAD](https://gnomad.broadinstitute.org), [CADD](https://cadd.gs.washington.edu) etc. ([`echtvar`](https://github.com/brentp/echtvar))
2. Annotate variants ([`VEP`](https://github.com/Ensembl/ensembl-vep))

##### Filtering

- TBD

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet with input data:

`samplesheet.csv`

```
sample,file,family_id,paternal_id,maternal_id,sex,phenotype
HG002,/path/to/HG002.fastq.gz,FAM1,HG003,HG004,1,1
HG005,/path/to/HG005.bam,FAM1,HG003,HG004,2,1
```

Now, you can run the pipeline using:

```bash
nextflow run genomic-medicine-sweden/skierfe -r dev -profile YOURPROFILE \
    --input samplesheet.csv \
    --preset <revio/pacbio/ONT_R10> \
    --outdir <OUTDIR> \
    --fasta <reference.fasta> \
    --skip_assembly_wf \
    --skip_repeat_wf \
    --skip_snv_annotation \
    --skip_cnv_calling
```

For more details and further functionality, please refer to the [usage documentation](https://github.com/genomic-medicine-sweden/skierfe/blob/dev/docs/usage.md).

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

To run in an offline environment, download the pipeline and singularity images using [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use):

```
nf-core download fellen31/skierfe -r dev
```

## Credits

genomic-medicine-sweden/skierfe was originally written by Felix Lenner.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use genomic-medicine-sweden/skierfe for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

> An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
