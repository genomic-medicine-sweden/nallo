# genomic-medicine-sweden/nallo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0dev - Unreleased

### `Added`

- [#148](https://github.com/genomic-medicine-sweden/nallo/pull/148) - Automatically infer sex if unknown
- [#148](https://github.com/genomic-medicine-sweden/nallo/pull/148) - Added read group tag to aligned BAM
- [#159](https://github.com/genomic-medicine-sweden/nallo/pull/159) - Allow files with from the same sample to be merged

### `Changed`

- [#146](https://github.com/genomic-medicine-sweden/nallo/pull/146) - Template merge for nf-core/tools v2.14.1
- [#145](https://github.com/genomic-medicine-sweden/nallo/pull/145) - Bump to new dev version
- [#151](https://github.com/genomic-medicine-sweden/nallo/pull/151) - Cleaned up TRGT output directory
- [#152](https://github.com/genomic-medicine-sweden/nallo/pull/152) - Use prefix in modkit module. Bgzip, index and split outputs into phased/unphased directories
- [#153](https://github.com/genomic-medicine-sweden/nallo/pull/153) - Changed cramino module to use prefix, renamed and moved all cramino outputs into `qc_alinged_reads/cramino/`
- [#159](https://github.com/genomic-medicine-sweden/nallo/pull/159) - Clarify the trio-binning genome assembly workflow
- [#159](https://github.com/genomic-medicine-sweden/nallo/pull/159) - `split_fastq` now splits on files instead of lines
- [#159](https://github.com/genomic-medicine-sweden/nallo/pull/159) - Use groupKey to remove bottleneck, where previously all samples had to wait before progressing after alignment

### `Fixed`

- [#156](https://github.com/genomic-medicine-sweden/nallo/pull/156) - Fixed program versions missing in output and MultiQC report

### Parameters

| Old parameter   | New parameter      |
| --------------- | ------------------ |
|                 | `--somalier_sites` |
| `--split_fastq` | `--split_fastq`\*  |

- `split_fastq` now splits the input files into _n_ files (range 2-999)

> [!NOTE]
> Parameter has been updated if both old and new parameter information is present.
> Parameter has been added if just the new parameter information is present.
> Parameter has been removed if new parameter information isn't present.

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
