# genomic-medicine-sweden/nallo: Output

## Aligned reads

[Minimap2](https://github.com/lh3/minimap2) is used to map the reads to a reference genome. The aligned reads are sorted, (merged) and indexed using [samtools](https://github.com/samtools/samtools).

| Path                                    | Description                         |
| --------------------------------------- | ----------------------------------- |
| `aligned_reads/minimap2/{sample}/*.bam` | Alignment file in bam format        |
| `aligned_reads/minimap2/{sample}/*.bai` | Index of the corresponding bam file |

If the pipeline is run with phasing, the aligned reads will be happlotagged using the active phasing tool.

| Path                                                              | Description             |
| ----------------------------------------------------------------- | ----------------------- |
| `{outputdir}/aligned_reads/{sample}/{sample}_haplotagged.bam`     | BAM file with haplotags |
| `{outputdir}/aligned_reads/{sample}/{sample}_haplotagged.bam.bai` | Index of the BAM file   |

!!!note

    Alignments will only be output without haplotags if phasing is off.

## Assembly

[Hifiasm](https://github.com/chhylp123/hifiasm) is used to assemble genomes. The assembled haplotypes are then comverted to fasta files using [gfastats](https://github.com/vgl-hub/gfastats). A deconstructed version of [dipcall](https://github.com/lh3/dipcall) is to map the assembled haplotypes back to the reference genome.

| Path                                                         | Description                                          |
| ------------------------------------------------------------ | ---------------------------------------------------- |
| `assembly_haplotypes/gfastats/{sample}/*hap1.p_ctg.fasta.gz` | Assembled haplotype 1                                |
| `assembly_haplotypes/gfastats/{sample}/*hap2.p_ctg.fasta.gz` | Assembled haplotype 2                                |
| `assembly_haplotypes/gfastats/{sample}/*.assembly_summary`   | Summary statistics                                   |
| `assembly_variant_calling/dipcall/{sample}/*hap1.bam`        | Assembled haplotype 1 mapped to the reference genome |
| `assembly_variant_calling/dipcall/{sample}/*hap1.bai`        | Index of the corresponding BAM file for haplotype 1  |
| `assembly_variant_calling/dipcall/{sample}/*hap2.bam`        | Assembled haplotype 2 mapped to the reference genome |
| `assembly_variant_calling/dipcall/{sample}/*hap2.bai`        | Index of the corresponding BAM file for haplotype 2  |

## Methylation pileups

[Modkit](https://github.com/nanoporetech/modkit) is used to create methylation pileups, producing bedMethyl files for both haplotagged and ungrouped reads. Additionaly, methylation information can be viewed in the BAM files, for example in IGV.

| Path                                                                                | Description                                               |
| ----------------------------------------------------------------------------------- | --------------------------------------------------------- |
| `methylation/modkit/pileup/phased/{sample}/*.modkit_pileup_phased_*.bed.gz`         | bedMethyl file with summary counts from haplotagged reads |
| `methylation/modkit/pileup/phased/{sample}/*.modkit_pileup_phased_ungrouped.bed.gz` | bedMethyl file for ungrouped reads                        |
| `methylation/modkit/pileup/unphased/{sample}/*.modkit_pileup.bed.gz`                | bedMethyl file with summary counts from all reads         |
| `methylation/modkit/pileup/unphased/{sample}/*.bed.gz.tbi`                          | Index of the corresponding bedMethyl file                 |

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
| `{outputdir}/aligned_reads/{sample}/{sample}_haplotagged.bam`     | BAM file with haplotags       |
| `{outputdir}/aligned_reads/{sample}/{sample}_haplotagged.bam.bai` | Index of the BAM file         |
| `{outputdir}/phased_variants/{sample}/*.vcf.gz`                   | VCF file with phased variants |
| `{outputdir}/phased_variants/{sample}/*.vcf.gz.tbi`               | Index of the VCF file         |
| `{outputdir}/qc/phasing_stats/{sample}/*.blocks.tsv`              | Phase block file              |
| `{outputdir}/qc/phasing_stats/{sample}/*.stats.tsv`               | Phasing statistics file       |

## QC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [cramino](https://github.com/wdecoster/cramino), [mosdepth](https://github.com/brentp/mosdepth), and [somalier](https://github.com/brentp/somalier) are used for read quality control.

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides general quality metrics for sequenced reads, including information on quality score distribution, per-base sequence content (%A/T/G/C), adapter contamination, and overrepresented sequences. For more details, refer to the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

| Path                                           | Description                                                     |
| ---------------------------------------------- | --------------------------------------------------------------- |
| `{outputdir}/qc/fastqc/{sample}/*_fastqc.html` | FastQC report containing quality metrics                        |
| `{outputdir}/qc/fastqc/{sample}/*_fastqc.zip`  | Zip archive with the FastQC report, data files, and plot images |

### Mosdepth

[Mosdepth](https://github.com/brentp/mosdepth) is used to report quality control metrics such as coverage and GC content from alignment files.

| Path                                                          | Description                                                                                                           |
| ------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `{outputdir}/qc/mosdepth/{sample}/*.mosdepth.global.dist.txt` | Cumulative distribution of bases covered for at least a given coverage value, across chromosomes and the whole genome |
| `{outputdir}/qc/mosdepth/{sample}/*.mosdepth.region.dist.txt` | Cumulative distribution of bases covered for at least a given coverage value, across regions (if a BED file is used)  |
| `{outputdir}/qc/mosdepth/{sample}/*.mosdepth.summary.txt`     | Mosdepth summary file                                                                                                 |
| `{outputdir}/qc/mosdepth/{sample}/*.regions.bed.gz`           | Depth per region (if a BED file is used)                                                                              |
| `{outputdir}/qc/mosdepth/{sample}/*.regions.bed.gz.csi`       | Index of the regions.bed.gz file                                                                                      |

### Cramino

[cramino](https://github.com/wdecoster/cramino) is used to analyze both phased and unphased reads.

| Path                                               | Description                                                                                          |
| -------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `{outputdir}/qc/cramino/phased/{sample}/*.arrow`   | Read length and quality in [Apache Arrow](https://arrow.apache.org/docs/format/Columnar.html) format |
| `{outputdir}/qc/cramino/phased/{sample}/*.txt`     | Summary information in text format                                                                   |
| `{outputdir}/qc/cramino/unphased/{sample}/*.arrow` | Read length and quality in [Apache Arrow](https://arrow.apache.org/docs/format/Columnar.html) format |
| `{outputdir}/qc/cramino/unphased/{sample}/*.txt`   | Summary information in text format                                                                   |

### Somalier

[somalier](https://github.com/brentp/somalier) checks relatedness and sex.

| Path                                                             | Description                                 |
| ---------------------------------------------------------------- | ------------------------------------------- |
| `{outputdir}/predigree/{project}.ped`                            | PED file updated with somalier-inferred sex |
| `{outputdir}/qc/somalier/relate/{project}/{project}.html`        | HTML report                                 |
| `{outputdir}/qc/somalier/relate/{project}/{project}.pairs.tsv`   | Information about sample pairs              |
| `{outputdir}/qc/somalier/relate/{project}/{project}.samples.tsv` | Information about individual samples        |

## Variants

### CNVs

[HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) is used to call CNVs. It also produces copy number, depth, and MAF [visualization tracks](#visualization-tracks).

| Path                                              | Description                               |
| ------------------------------------------------- | ----------------------------------------- |
| `svs/family/{family_id}/${family_id_cnvs_merged..vcf.gz`           | VCF file containing CNV variants          |
| `svs/family/{family_id}/{family_id}_cnvs_merged.vcf.gz.tbi`       | Index of the corresponding VCF file       |

### Paralogous genes

[Paraphase](https://github.com/PacificBiosciences/paraphase) is used to call paralogous genes.

| Path                                                        | Description                             |
| ----------------------------------------------------------- | --------------------------------------- |
| `paraphase/{sample}/*.bam`                                  | BAM file with haplotypes grouped by HP  |
| `paraphase/{sample}/*.bai`                                  | Index of the BAM file                   |
| `paraphase/{sample}/*.json`                                 | Summary of haplotypes and variant calls |
| `paraphase/{sample}_paraphase_vcfs/{sample}_{gene}_vcf`     | VCF file per gene                       |
| `paraphase/{sample}_paraphase_vcfs/{sample}_{gene}_vcf.tbi` | Index of the VCF file                   |

### Repeats

[TRGT](https://github.com/PacificBiosciences/trgt) is used to call repeats:

| Path                                                                  | Description                               |
| --------------------------------------------------------------------- | ----------------------------------------- |
| `{outputdir}/repeat_calling/trgt/multi_sample/{project}/*.vcf.gz`     | Merged VCF file for all samples           |
| `{outputdir}/repeat_calling/trgt/multi_sample/{project}/*.vcf.gz.tbi` | Index of the VCF file                     |
| `{outputdir}/repeat_calling/trgt/single_sample/{sample}/*.vcf.gz`     | VCF file with called repeats for a sample |
| `{outputdir}/repeat_calling/trgt/single_sample/{sample}/*.vcf.gz.tbi` | Index of the VCF file                     |
| `{outputdir}/repeat_calling/trgt/single_sample/{sample}/*.bam`        | BAM file with sorted spanning reads       |
| `{outputdir}/repeat_calling/trgt/single_sample/{sample}/*.bai`        | Index of the BAM file                     |

[Stranger](https://github.com/Clinical-Genomics/stranger) is used to annotate them:

| Path                                                           | Description                     |
| -------------------------------------------------------------- | ------------------------------- |
| `{outputdir}/repeat_annotation/stranger/{sample}/*.vcf.gz`     | Annotated VCF file              |
| `{outputdir}/repeat_annotation/stranger/{sample}/*.vcf.gz.tbi` | Index of the annotated VCF file |

### SNVs

[DeepVariant](https://github.com/google/deepvariant) is used to call variants, while [bcftools](https://samtools.github.io/bcftools/bcftools.html) and [GLnexus](https://github.com/dnanexus-rnd/GLnexus) are used for merging variants.

!!!note

    Variants are only output without annotation and ranking if these subworkflows are turned off.

| Path                                                   | Description                                                                 |
| ------------------------------------------------------ | --------------------------------------------------------------------------- |
| `snvs/single_sample/{sample}/{sample}_snv.vcf.gz`      | VCF file containing called variants with alternative genotypes for a sample |
| `snvs/single_sample/{sample}/{sample}_snv.vcf.gz.tbi`  | Index of the corresponding VCF file                                         |
| `snvs/multi_sample/{project}/{project}_snv.vcf.gz`     | VCF file containing called variants for all samples                         |
| `snvs/multi_sample/{project}/{project}_snv.vcf.gz.tbi` | Index of the corresponding VCF file                                         |
| `snvs/stats/single_sample/*.stats.txt`                 | Variant statistics                                                          |

[echtvar](https://github.com/brentp/echtvar) and [VEP](https://www.ensembl.org/vep) are used for annotating SNVs, while [CADD](https://cadd.gs.washington.edu/) is used to annotate INDELs with CADD scores.

!!!note

    Variants are only output without ranking if that subworkflows are turned off.

| Path                                                             | Description                                                                    |
| ---------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| `databases/echtvar/encode/{project}/*.zip`                       | Database with allele frequency (AF) and allele count (AC) for all samples      |
| `snvs/single_sample/{sample}/{sample}_snv_annotated.vcf.gz`      | VCF file containing annotated variants with alternative genotypes for a sample |
| `snvs/single_sample/{sample}/{sample}_snv_annotated.vcf.gz.tbi`  | Index of the annotated VCF file                                                |
| `snvs/multi_sample/{project}/{project}_snv_annotated.vcf.gz`     | VCF file containing annotated variants for all samples                         |
| `snvs/multi_sample/{project}/{project}_snv_annotated.vcf.gz.tbi` | Index of the annotated VCF file                                                |

[GENMOD](https://github.com/Clinical-Genomics/genmod) is used to rank the annotated SNVs and INDELs.

| Path                                                                    | Description                                                 |
| ----------------------------------------------------------------------- | ----------------------------------------------------------- |
| `snvs/single_sample/{sample}/{sample}_snv_annotated_ranked.vcf.gz`      | VCF file with annotated and ranked variants for a sample    |
| `snvs/single_sample/{sample}/{sample}_snv_annotated_ranked.vcf.gz.tbi`  | Index of the ranked VCF file                                |
| `snvs/multi_sample/{project}/{project}_snv_annotated_ranked.vcf.gz`     | VCF file with annotated and ranked variants for all samples |
| `snvs/multi_sample/{project}/{project}_snv_annotated_ranked.vcf.gz.tbi` | Index of the ranked VCF file                                |

### SVs

[Severus](https://github.com/KolmogorovLab/Severus) or [Sniffles](https://github.com/fritzsedlazeck/Sniffles) is used to call structural variants, and [SVDB](https://github.com/J35P312/SVDB) is used to merge variants within and between samples.

!!!note

    Variants are only output without annotation if that subworkflow is turned off.

| Path                                                  | Description                                                  |
| ----------------------------------------------------- | ------------------------------------------------------------ |
| `svs/multi_sample/{project}/{project}_svs.vcf.gz`     | VCF file with merged structural variants for all samples     |
| `svs/multi_sample/{project}/{project}_svs.vcf.gz.tbi` | Index of the merged VCF file                                 |
| `svs/single_sample/{sample}/*.vcf.gz`                 | VCF file with merged structural variants for a single sample |
| `svs/single_sample/{sample}/*.vcf.gz.tbi`             | Index of the VCF file                                        |

[SVDB](https://github.com/J35P312/SVDB) and [VEP](https://www.ensembl.org/vep) are used to annotate structural variants.

| Path                                                            | Description                                                        |
| --------------------------------------------------------------- | ------------------------------------------------------------------ |
| `svs/multi_sample/{project}/{project}_svs_annotated.vcf.gz`     | VCF file with annotated merged structural variants for all samples |
| `svs/multi_sample/{project}/{project}_svs_annotated.vcf.gz.tbi` | Index of the annotated VCF file                                    |
| `svs/single_sample/{sample}/*.vcf_annotated.gz`                 | VCF file with annotated structural variants for a single sample    |
| `svs/single_sample/{sample}/*.vcf_annotated.gz.tbi`             | Index of the annotated VCF file                                    |

## Visualization Tracks

[HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) is used to call CNVs, but it also produces copy number, depth, and MAF tracks that can be visualized in for example IGV.

| `visualization_tracks/{sample}/*.copynum.bedgraph` | Copy number in bedgraph format            |
| `visualization_tracks/{sample}/*.depth.bw`         | Depth track in BigWig format              |
| `visualization_tracks/{sample}/*.maf.bw`           | Minor allele frequencies in BigWig format |
