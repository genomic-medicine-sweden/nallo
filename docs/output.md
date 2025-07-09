# genomic-medicine-sweden/nallo: Output

## Introduction

This document describes the pipeline output files and the tools used to generate them.

## Aligned reads

[Minimap2](https://github.com/lh3/minimap2) is used to map the reads to a reference genome. The aligned reads are sorted, merged and indexed using [samtools](https://github.com/samtools/samtools). If the pipeline is run with phasing, the aligned reads will be haplotagged using the active phasing tool.

| Path                                                                   | Description                               | Alignment          | Alignment & phasing |
| ---------------------------------------------------------------------- | ----------------------------------------- | ------------------ | ------------------- |
| `aligned_reads/minimap2/{sample}/{sample}_aligned.{bam,cram}`          | Alignment file in BAM or CRAM format      | :white_check_mark: |                     |
| `aligned_reads/minimap2/{sample}/{sample}_aligned.{bam.bai,cram.crai}` | Index of the corresponding alignment file | :white_check_mark: |                     |

| Path                                                              | Description                                         | Alignment | Alignment & phasing |
| ----------------------------------------------------------------- | --------------------------------------------------- | --------- | ------------------- |
| `aligned_reads/{sample}/{sample}_haplotagged.{bam,cram}`          | Alignment file with haplotags in BAM or CRAM format |           | :white_check_mark:  |
| `aligned_reads/{sample}/{sample}_haplotagged.{bam.bai,cram.crai}` | Index of the alignment file                         |           | :white_check_mark:  |

## Assembly

[Hifiasm](https://github.com/chhylp123/hifiasm) is used to assemble genomes. The assembled haplotypes are then aligned to the reference genome with [minimap2](https://github.com/lh3/minimap2), tagged with `HP:1` for the "paternal" haplotype, and `HP:2` for the "maternal" haplotype, before being merged together into one file with [samtools](https://github.com/samtools/samtools). [gfastats](https://github.com/vgl-hub/gfastats) is used to convert the assembly to fasta format before alignment, and also outputs summary stats per haplotype.

| Path                                                                     | Description                                                                                       |
| ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------- |
| `assembly/sample/{sample}/{sample}_aligned_assembly.{bam,cram}`          | Both assembled haplotypes mapped to the reference genome, merged and haplotagged (`HP:1`/`HP:2`). |
| `assembly/sample/{sample}/{sample}_aligned_assembly.{bam.bai,cram.crai}` | Index of aligned assembly.                                                                        |
| `assembly/stats/{sample}/{sample}_haplotype_1.assembly_summary`          | Summary statistics for haplotype 1/paternal haplotype                                             |
| `assembly/stats/${sample}/{sample}_haplotype_2.assembly_summary`         | Summary statistics for haplotype 2/maternal haplotype                                             |

## Methylation pileups

[Modkit](https://github.com/nanoporetech/modkit) is used to create methylation pileups, producing bedMethyl files for both haplotagged and ungrouped reads. Additionally, methylation information can be viewed in the BAM files, for example in IGV. When phasing is on, modkit outputs pileups per haplotype.

| Path                                                                  | Description                                                             | Alignment          | Alignment & phasing |
| --------------------------------------------------------------------- | ----------------------------------------------------------------------- | ------------------ | ------------------- |
| `methylation/modkit/pileup/{sample}/*.modkit_pileup_1.bed.gz`         | bedMethyl file with summary counts from haplotagged reads (haplotype 1) |                    | :white_check_mark:  |
| `methylation/modkit/pileup/{sample}/*.modkit_pileup_2.bed.gz`         | bedMethyl file with summary counts from haplotagged reads (haplotype 2) |                    | :white_check_mark:  |
| `methylation/modkit/pileup/{sample}/*.modkit_pileup_ungrouped.bed.gz` | bedMethyl file for ungrouped reads                                      |                    | :white_check_mark:  |
| `methylation/modkit/pileup/{sample}/*.modkit_pileup.bed.gz`           | bedMethyl file with summary counts from all reads                       | :white_check_mark: |                     |
| `methylation/modkit/pileup/{sample}/*.bed.gz.tbi`                     | Index of the corresponding bedMethyl files                              | :white_check_mark: |                     |

## MultiQC

[MultiQC](http://multiqc.info) generates an HTML report summarizing all samples' QC results and pipeline statistics.

| Path                          | Description                               |
| ----------------------------- | ----------------------------------------- |
| `multiqc/multiqc_report.html` | HTML report summarizing QC results        |
| `multiqc/multiqc_data/`       | Directory containing parsed statistics    |
| `multiqc/multiqc_plots/`      | Directory containing static report images |

## Pipeline Information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) generates reports for troubleshooting, performance, and traceability.

| Path                                    | Description                       |
| --------------------------------------- | --------------------------------- |
| `pipeline_info/execution_report.html`   | Execution report                  |
| `pipeline_info/execution_timeline.html` | Timeline report                   |
| `pipeline_info/execution_trace.txt`     | Execution trace                   |
| `pipeline_info/pipeline_dag.dot`        | Pipeline DAG in DOT format        |
| `pipeline_info/pipeline_report.html`    | Pipeline report                   |
| `pipeline_info/software_versions.yml`   | Software versions used in the run |

## Phasing

[LongPhase](https://github.com/twolinin/longphase), [WhatsHap](https://whatshap.readthedocs.io/en/latest/), or [HiPhase](https://github.com/PacificBiosciences/HiPhase) are used for phasing.

| Path                                                              | Description                   |
| ----------------------------------------------------------------- | ----------------------------- |
| `aligned_reads/{sample}/{sample}_haplotagged.{bam,cram}`          | BAM/CRAM file with haplotags  |
| `aligned_reads/{sample}/{sample}_haplotagged.{bam.bai,cram.crai}` | Index of the BAM/CRAM file    |
| `phased_variants/{sample}/*.vcf.gz`                               | VCF file with phased variants |
| `phased_variants/{sample}/*.vcf.gz.tbi`                           | Index of the VCF file         |
| `qc/phasing_stats/{sample}/*.blocks.gtf.gz`                       | Phase block file              |
| `qc/phasing_stats/{sample}/*.blocks.gtf.gz.tbi`                   | Index of block file           |
| `qc/phasing_stats/{sample}/*.stats.tsv`                           | Phasing statistics file       |

## QC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [cramino](https://github.com/wdecoster/cramino), [mosdepth](https://github.com/brentp/mosdepth), and [somalier](https://github.com/brentp/somalier) are used for read quality control.

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides general quality metrics for sequenced reads, including information on quality score distribution, per-base sequence content (%A/T/G/C), adapter contamination, and overrepresented sequences. For more details, refer to the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

| Path                               | Description                                                     |
| ---------------------------------- | --------------------------------------------------------------- |
| `qc/fastqc/{sample}/*_fastqc.html` | FastQC report containing quality metrics                        |
| `qc/fastqc/{sample}/*_fastqc.zip`  | Zip archive with the FastQC report, data files, and plot images |

### Mosdepth

[Mosdepth](https://github.com/brentp/mosdepth) is used to report quality control metrics such as coverage and GC content from alignment files.

| Path                                                     | Description                                                                                                           | With `--target_regions` | Without `--target_regions` |
| -------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- | ----------------------- | -------------------------- |
| `qc/mosdepth/{sample}/{sample}.mosdepth.global.dist.txt` | Cumulative distribution of bases covered for at least a given coverage value, across chromosomes and the whole genome | :white_check_mark:      | :white_check_mark:         |
| `qc/mosdepth/{sample}/{sample}.mosdepth.summary.txt`     | Mosdepth summary file                                                                                                 | :white_check_mark:      | :white_check_mark:         |
| `qc/mosdepth/{sample}/{sample}.mosdepth.region.dist.txt` | Cumulative distribution of bases covered for at least a given coverage value, across regions                          | :white_check_mark:      |                            |
| `qc/mosdepth/{sample}/{sample}.per-base.d4`              | Per-base depth in d4 format                                                                                           | :white_check_mark:      |
| `qc/mosdepth/{sample}/{sample}.regions.bed.gz`           | Depth per region                                                                                                      | :white_check_mark:      |
| `qc/mosdepth/{sample}/{sample}.regions.bed.gz.csi`       | Index of the regions.bed.gz file                                                                                      | :white_check_mark:      |

### Cramino

[cramino](https://github.com/wdecoster/cramino) is used to analyze both phased and unphased reads.

| Path                                   | Description                                                                                          |
| -------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `qc/cramino/phased/{sample}/*.arrow`   | Read length and quality in [Apache Arrow](https://arrow.apache.org/docs/format/Columnar.html) format |
| `qc/cramino/phased/{sample}/*.txt`     | Summary information in text format                                                                   |
| `qc/cramino/unphased/{sample}/*.arrow` | Read length and quality in [Apache Arrow](https://arrow.apache.org/docs/format/Columnar.html) format |
| `qc/cramino/unphased/{sample}/*.txt`   | Summary information in text format                                                                   |

### Somalier

[somalier](https://github.com/brentp/somalier) checks relatedness and sex.

| Path                                                 | Description                                            |
| ---------------------------------------------------- | ------------------------------------------------------ |
| `pedigree/family/{family).ped`                       | PED file updated with somalier-inferred sex per family |
| `qc/somalier/relate/{project}/{project}.html`        | HTML report                                            |
| `qc/somalier/relate/{project}/{project}.pairs.tsv`   | Information about sample pairs                         |
| `qc/somalier/relate/{project}/{project}.samples.tsv` | Information about individual samples                   |

### Peddy

[peddy](https://github.com/brentp/peddy) checks relatedness and sex.

| Path                                                      | Description                                                                                                                                                                   |
| --------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `qc/peddy/{family}/{family}.peddy.ped`                    | PED file updated with peddy-inferred sex per family                                                                                                                           |
| `qc/peddy/{family}/{family}.html`                         | HTML report                                                                                                                                                                   |
| `qc/peddy/{family}/{family}.vs.html`                      | HTML report of observed vs expected relatedness                                                                                                                               |
| `qc/peddy/{family}/{family}.sex_check.csv`                | Comparison between reported sex (ped file) and that inferred from peddy                                                                                                       |
| `qc/peddy/{family}/{family}.het_check.csv`                | Het check does general QC including rate of het calls, allele-balance at het calls, mean and median depth, and a PCA projection onto thousand genomes. Incudes ancestry check |
| `qc/peddy/{family}/{family}.ped_check.csv`                | Ped check compares the relatedness of 2 samples as reported in a .ped file to the relatedness inferred from the genotypes and ~25K sites in the genome                        |
| `qc/peddy/{family}/{family}.sex_check.png`                | PNG comparison between reported sex (ped file) and that inferred from peddy                                                                                                   |
| `qc/peddy/{family}/{family}.het_check.png`                | PNG of heterozygosity check                                                                                                                                                   |
| `qc/peddy/{family}/{family}.ped_check.png`                | PNG of the ped check comparison                                                                                                                                               |
| `qc/peddy/{family}/{family}.ped_check.rel-difference.csv` | CSV file with the comparison between inferred and given relatedness                                                                                                           |

### DeepVariant

`vcf_stats_report.py` from [DeepVariant](https://github.com/google/deepvariant) is used to generate a html report per sample.

| Path                                                                  | Description                                 |
| --------------------------------------------------------------------- | ------------------------------------------- |
| `qc/deepvariant_vcfstatsreport/{sample}/${sample}.visual_report.html` | Visual report of SNV calls from DeepVariant |

### Bcftools

[bcftools stats](https://samtools.github.io/bcftools/bcftools.html) is used to generate variants statistics from SNV calls.

| Path                                                       | Description        |
| ---------------------------------------------------------- | ------------------ |
| `qc/bcftools_stats/${sample}/${sample}.bcftools_stats.txt` | Variant statistics |

## Variants

In general, annotated variant calls are output per family while unannotated calls are output per sample.

### Paralogous genes

[Paraphase](https://github.com/PacificBiosciences/paraphase) is used to call paralogous genes.

| Path                                                                   | Description                                               |
| ---------------------------------------------------------------------- | --------------------------------------------------------- |
| `paraphase/sample/{sample}/{sample}.paraphase.{bam,cram}`              | BAM/CRAM file with reads from analyzed regions            |
| `paraphase/sample/{sample}/{sample}.paraphase.{bam.bai,cram.crai}`     | Index of the BAM/CRAM file                                |
| `paraphase/sample/{sample}/{sample}.paraphase.json`                    | Summary of haplotypes and variant calls                   |
| `paraphase/sample/{sample}_paraphase_vcfs/{sample}_{gene}_vcf.gz`      | VCF file per gene                                         |
| `paraphase/sample/{sample}_paraphase_vcfs/{sample}_{gene}_vcf.gz.tbi`  | Index of the VCF file                                     |
| `paraphase/family/{family_id}/{family_id}_paraphase_merged.vcf.gz`     | VCF file from paraphase, merged by family                 |
| `paraphase/family/{family_id}/{family_id}_paraphase_merged.vcf.gz.tbi` | Index of the VCF file merged by family                    |
| `paraphase/family/{family_id}/{family_id}_merged.json`                 | Summary of haplotypes and variant calls, merged by family |

### Repeats

[TRGT](https://github.com/PacificBiosciences/trgt) or [STRdust](https://github.com/wdecoster/STRdust) are used to call repeats.
[Stranger](https://github.com/Clinical-Genomics/stranger) is used to annotate repeats.

| Path                                                                      | Description                               | STRdust, Call repeats only | TRGT, Call repeats only | TRGT, Call & annotate repeats |
| ------------------------------------------------------------------------- | ----------------------------------------- | -------------------------- | ----------------------- | ----------------------------- |
| `repeats/sample/{sample}/{sample}_{str_caller}.vcf.gz`                    | VCF file with called repeats for a sample | :white_check_mark:         | :white_check_mark:      | :white_check_mark:            |
| `repeats/sample/{sample}/{sample}_{str_caller}.vcf.gz.tbi`                | Index of the VCF file                     | :white_check_mark:         | :white_check_mark:      | :white_check_mark:            |
| `repeats/sample/{sample}/{sample}_spanning_trgt.{bam,cram}`               | BAM/CRAM file with sorted spanning reads  |                            | :white_check_mark:      | :white_check_mark:            |
| `repeats/sample/{sample}/{sample}_spanning_trgt.{bam.bai,cram.crai}`      | Index of the BAM/CRAM file                |                            | :white_check_mark:      | :white_check_mark:            |
| `repeats/family/{family}/{family}_repeat_expansions.vcf.gz`               | Merged VCF file per family                | :white_check_mark:         | :white_check_mark:      |                               |
| `repeats/family/{family}/{family}_repeat_expansions.vcf.gz.tbi`           | Index of the VCF file                     | :white_check_mark:         | :white_check_mark:      |                               |
| `repeats/family/{family}/{family}_repeat_expansions_annotated.vcf.gz`     | Merged, annotated VCF file per family     |                            |                         | :white_check_mark:            |
| `repeats/family/{family}/{family}_repeat_expansions_annotated.vcf.gz.tbi` | Index of the VCF file                     |                            |                         | :white_check_mark:            |

### SNVs

[DeepVariant](https://github.com/google/deepvariant) is used to call variants, while [bcftools](https://samtools.github.io/bcftools/bcftools.html) and [GLnexus](https://github.com/dnanexus-rnd/GLnexus) are used for merging variants.

| Path                                                                  | Description                                                     | Call SNVs          | Call & annotate SNVs | Call, annotate and rank SNVs |
| --------------------------------------------------------------------- | --------------------------------------------------------------- | ------------------ | -------------------- | ---------------------------- |
| `qc/bcftools_stats/${sample}/${sample}.bcftools_stats.txt`            | Variant statistics                                              | :white_check_mark: | :white_check_mark:   | :white_check_mark:           |
| `qc/deepvariant_vcfstatsreport/{sample}/${sample}.visual_report.html` | Visual report of SNV calls from DeepVariant                     | :white_check_mark: | :white_check_mark:   | :white_check_mark:           |
| `snvs/sample/{sample}/{sample}_deepvariant_snvs.vcf.gz`               | VCF file containing called variants for a single sample         | :white_check_mark: | :white_check_mark:   | :white_check_mark:           |
| `snvs/sample/{sample}/{sample}_deepvariant_snvs.vcf.gz.tbi`           | Index of the corresponding VCF file                             | :white_check_mark: | :white_check_mark:   | :white_check_mark:           |
| `snvs/family/{family}/{family}_snvs.vcf.gz`                           | VCF file containing called variants for all samples in a family | :white_check_mark: |                      |                              |
| `snvs/family/{family}/{family}_snvs.vcf.gz.tbi`                       | Index of the corresponding VCF file                             | :white_check_mark: |                      |                              |

#### Annotation

[Echtvar](https://github.com/brentp/echtvar) and [VEP](https://www.ensembl.org/vep) are used for annotating SNVs, while [CADD](https://cadd.gs.washington.edu/) is used to annotate INDELs with CADD scores.

| Path                                                               | Description                                       | Call SNVs | Call & annotate SNVs | Call, annotate and rank SNVs |
| ------------------------------------------------------------------ | ------------------------------------------------- | --------- | -------------------- | ---------------------------- |
| `snvs/family/{family}/{family}_snvs_annotated_research.vcf.gz`     | VCF file containing annotated variants per family |           | :white_check_mark:   |                              |
| `snvs/family/{family}/{family}_snvs_annotated_research.vcf.gz.tbi` | Index of the annotated VCF file                   |           | :white_check_mark:   |                              |

#### Ranking

[GENMOD](https://github.com/Clinical-Genomics/genmod) is used to rank the annotated SNVs and INDELs.

| Path                                                                      | Description                                            | Call SNVs | Call & annotate SNVs | Call, annotate and rank SNVs |
| ------------------------------------------------------------------------- | ------------------------------------------------------ | --------- | -------------------- | ---------------------------- |
| `snvs/family/{family}/{family}_snvs_annotated_ranked_research.vcf.gz`     | VCF file with annotated and ranked variants per family |           |                      | :white_check_mark:           |
| `snvs/family/{family}/{family}_snvs_annotated_ranked_research.vcf.gz.tbi` | Index of the ranked VCF file                           |           |                      | :white_check_mark:           |

#### Filtering

[Filter_vep](https://www.ensembl.org/vep) and [bcftools](https://samtools.github.io/bcftools/bcftools.html) can be used to filter variants after annotation. These will be output if either of `--filter_variants_hgnc_id` and `--filter_snvs_expression` has been used, and only family VCFs are filtered.

| Path                                           | Description                                  |
| ---------------------------------------------- | -------------------------------------------- |
| `snvs/{family}/{family}_*_clinical.vcf.gz`     | VCF file with filtered variants for a family |
| `snvs/{family}/{family}_*_clinical.vcf.gz.tbi` | Index of the filtered VCF file               |

!!!tip

    Filtered variants are output alongside unfiltered variants as additional files.

### SVs

[Severus](https://github.com/KolmogorovLab/Severus) or [Sniffles](https://github.com/fritzsedlazeck/Sniffles) are used to call structural variants, while [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) is used to call CNVs. HiFiCNV also produces copy number, depth, and MAF [visualization tracks](#visualization-tracks).

!!!tip "Family-level VCFs per caller"

    Unannotated family-level VCFs per caller can be output with --publish_unannotated_family_svs.

| Path                                                                           | Description                                        | Call SVs           | `--publish_unannotated_family_svs` |
| ------------------------------------------------------------------------------ | -------------------------------------------------- | ------------------ | ---------------------------------- | --- |
| `svs/sample/{family_id}/{family_id}_{sniffles,severus,hificnv}_svs.vcf.gz`     | VCF file with merged SVs/CNVs by family and caller |                    | :white_check_mark:                 |     |
| `svs/sample/{family_id}/{family_id}_{sniffles,severus,hificnv}_svs.vcf.gz.tbi` | Index of the merged VCF file                       |                    | :white_check_mark:                 |     |
| `svs/family/{family_id}/{family_id}_svs.vcf.gz`                                | VCF file with merged SVs/CNVs by family            | :white_check_mark: |                                    |
| `svs/family/{family_id}/{family_id}_svs.vcf.gz.tbi`                            | Index of the merged VCF file                       | :white_check_mark: |                                    |

#### Annotation

[SVDB](https://github.com/J35P312/SVDB) and [VEP](https://www.ensembl.org/vep) are used to annotate structural variants.

| Path                                                                   | Description                                                |
| ---------------------------------------------------------------------- | ---------------------------------------------------------- |
| `svs/family/{family_id}/{family_id}_svs_annotated_research.vcf.gz`     | VCF file with merged and annotated CNVs and SVs per family |
| `svs/family/{family_id}/{family_id}_svs_annotated_research.vcf.gz.tbi` | Index of the merged VCF file                               |

#### Ranking

[GENMOD](https://github.com/Clinical-Genomics/genmod) is used to rank the annotated SVs.

| Path                                                                          | Description                                                        |
| ----------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| `svs/family/{family_id}/{family_id}_svs_annotated_ranked_research.vcf.gz`     | VCF file with merged, annotated and ranked CNVs and SVs per family |
| `svs/family/{family_id}/{family_id}_svs_annotated_ranked_research.vcf.gz.tbi` | Index of the merged VCF file                                       |

#### Filtering

[Filter_vep](https://www.ensembl.org/vep) and [bcftools](https://samtools.github.io/bcftools/bcftools.html) can be used to filter variants after annotation. These will be output if either of `--filter_variants_hgnc_id` and `--filter_svs_expression` has been used, and only family VCFs are filtered.

| Path                                          | Description                                  |
| --------------------------------------------- | -------------------------------------------- |
| `svs/{family}/{family}_*_clinical.vcf.gz`     | VCF file with filtered variants for a family |
| `svs/{family}/{family}_*_clinical.vcf.gz.tbi` | Index of the filtered VCF file               |

!!!tip

    Filtered variants are output alongside unfiltered variants as additional files.

## Visualization Tracks

[HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) is used to call CNVs, but it also produces copy number, depth, and MAF tracks that can be visualized in for example IGV.

| Path                                               | Description                               |
| -------------------------------------------------- | ----------------------------------------- |
| `visualization_tracks/{sample}/*.copynum.bedgraph` | Copy number in bedgraph format            |
| `visualization_tracks/{sample}/*.depth.bw`         | Depth track in BigWig format              |
| `visualization_tracks/{sample}/*.maf.bw`           | Minor allele frequencies in BigWig format |
