# genomic-medicine-sweden/nallo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0dev - Unreleased

### `Added`

### `Changed`

- Template merge for nf-core/tools v2.14.1 [#146](https://github.com/genomic-medicine-sweden/nallo/pull/146)
- Bump to new dev version [#145](https://github.com/genomic-medicine-sweden/nallo/pull/145)

### `Fixed`

## v0.1.0 - [2024-05-08]

Initial release of genomic-medicine-sweden/nallo, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Raw read QC with FastQC and FQCRS
- Align reads to reference with minimap2
- Aligned read QC with cramino and mosdepth
- Call SNVs with DeepVariant and merge with GLNexus
- Annotate SNVs with echtvar and VEP
- Call SVs with Sniffles, tandem repeats with TRGT and CNVs with HiFiCNV
- Phase variants and haplotag reads with whatshap or HiPhase
- Create methylation pileups with modkit
- Assemble genomes with hifiasm
- Align assemly to reference and call variants with dipcall
