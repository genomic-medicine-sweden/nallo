

# genomic-medicine-sweden/nallo pipeline parameters                                                                                
                                                                                                                                   
Long-read variant calling pipeline                                                                                                 
                                                                                                                                   
## Workflow skip options                                                                                                           
                                                                                                                                   
Allows skipping certain parts of the pipeline                                                                                      
                                                                                                                                   
| Parameter | Description | Type | Default | Required | Hidden |                                                                   
|-----------|-----------|-----------|-----------|-----------|-----------|                                                          
| `skip_qc` | Skip QC of reads | `boolean` | False |  |  |                                                                         
| `skip_snv_calling` | Skip short variant calling | `boolean` | False |  |  |                                                      
| `skip_genome_assembly` | Skip genome assembly and assembly variant calling | `boolean` | False |  |  |                           
| `skip_alignment` | Skip read mapping (alignment) | `boolean` | False |  |  |                                                     
| `skip_methylation_calling` | Skip methylation calling | `boolean` | False |  |  |                                                
| `skip_repeat_calling` | Skip tandem repeat calling | `boolean` | False |  |  |                                                   
| `skip_repeat_annotation` | Skip tandem repeat annotation | `boolean` | False |  |  |                                             
| `skip_chromograph` | Skip chromograph image generation. True if both plot_chromograph_coverage and plot_chromograph_autozygosity 
are set to false. | `boolean` |  |  |  |                                                                                           
| `skip_peddy` | Skip peddy | `boolean` | False |  |  |                                                                            
| `skip_sambamba_depth` | Skip sambamba depth. True unless `--sambamba_regions` is provided. | `boolean` |  |  |  |                
| `skip_phasing` | Skip phasing of variants and haplotagging of reads | `boolean` | False |  |  |                                  
| `skip_snv_annotation` | Skip short variant annotation | `boolean` | False |  |  |                                                
| `skip_sv_calling` | Skip structural variant calling | `boolean` | False |  |  |                                                  
| `skip_sv_annotation` | Skip structural variant annotation | `boolean` | False |  |  |                                            
| `skip_call_paralogs` | Skip the calling of specific paralogous genes | `boolean` | False |  |  |                                 
| `skip_rank_variants` | Skip ranking of short variants | `boolean` | False |  |  |                                                
| `skip_prepare_gens_input` | Skip preparing input data for Gens | `boolean` | False |  |  |                                       
                                                                                                                                   
## Input/output options                                                                                                            
                                                                                                                                   
Define where the pipeline should find input data and save output data.                                                             
                                                                                                                                   
| Parameter | Description | Type | Default | Required | Hidden |                                                                   
|-----------|-----------|-----------|-----------|-----------|-----------|                                                          
| `input` | Path to comma-separated file containing information about the samples in the experiment.                               
<details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment 
before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a 
header row.</small></details>| `string` |  | True |  |                                                                             
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud              
infrastructure. | `string` |  | True |  |                                                                                          
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address 
to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file               
(`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  
|                                                                                                                                  
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  | 
|                                                                                                                                  
| `cadd_prescored_indels` | Path to a directory containing prescored indels for CADD. <details><summary>Help</summary><small>This  
folder contains the compressed files and indexes that would otherwise be in data/prescored folder as described in                  
https://github.com/kircherlab/CADD-scripts/#manual-installation.</small></details>| `string` |  |  |  |                            
| `cadd_resources` | Path to a directory containing CADD annotations. <details><summary>Help</summary><small>This folder contains  
the uncompressed files that would otherwise be in data/annotation folder as described in                                           
https://github.com/kircherlab/CADD-scripts/#manual-installation.</small></details>| `string` |  |  |  |                            
| `par_regions` | Provide a bed file of chrX and chrY PAR regions for DeepVariant | `string` |  |  |  |                            
| `tandem_repeats` | A tandem repeat BED file for sniffles | `string` |  |  |  |                                                   
| `str_bed` | A BED file with repeats to be genotyped with TRGT | `string` |  |  |  |                                              
| `echtvar_snv_databases` | Path to a CSV/TSV/JSON/YAML file with echtvar databases to annotate SNVs with. | `string` |  |  |  |   
| `svdb_sv_databases` | Databases used for structural variant annotation in vcf format. <details><summary>Help</summary><small>Path
to a CSV/TSV/JSON/YAML file containing information about the databases used for structural variant annotation.</small></details>|  
`string` |  |  |  |                                                                                                                
| `stranger_repeat_catalog` | A variant catalog json-file for stranger | `string` |  |  |  |                                       
| `variant_consequences_snvs` | File containing list of SO terms listed in the order of severity from most severe to lease severe  
for annotating genomic SNVs. For more information check https://ensembl.org/info/genome/variation/prediction/predicted_data.html | 
`string` |  |  |  |                                                                                                                
| `variant_consequences_svs` | File containing list of SO terms listed in the order of severity from most severe to lease severe   
for annotating genomic SVs. For more information check https://ensembl.org/info/genome/variation/prediction/predicted_data.html |  
`string` |  |  |  |                                                                                                                
| `vep_cache` | A path to the VEP cache location | `string` |  |  |  |                                                             
| `target_regions` | A BED file with regions of interest. | `string` |  |  |  |                                                    
| `methbat_regions` | A tsv file with only regions of interest (example here:                                                      
https://github.com/PacificBiosciences/MethBat/blob/main/data/cpgIslandExt.sorted.hg38.tsv) or with both regions and background     
cohort values (example here: https://github.com/PacificBiosciences/MethBat/blob/main/data/meth_profile_model.tsv), made with       
methbat build | `string` |  |  |  |                                                                                                
| `modkit_call_regions` | A BED file with regions of interest for the methylation pileups. By default this is the same as          
`target_regions`. | `string` |  |  |  |                                                                                            
| `mosdepth_regions` | A BED file with regions of interest used in mosdepth. By default this is the same as `qc_regions`. |        
`string` |  |  |  |                                                                                                                
| `mosdepth_d4_output` | Should mosdepth output d4-files. | `boolean` | False |  |  |                                              
| `bigwig_modcodes` | Comma-separated list of modification codes to include in the bigWig methylation visualization file. Defaults 
to 5hmC and 5mC. See https://samtools.github.io/hts-specs/SAMtags.pdf for a complete list. | `string` | h,m |  |  |                
| `snv_call_regions` | A BED file containing regions to limit SNV calling. By default this is the same as `target_regions`. |      
`string` |  |  |  |                                                                                                                
| `sv_call_regions` | A BED file containging regions to filter SV calls. By default this is the same as `target_regions`. |        
`string` |  |  |  |                                                                                                                
| `qc_regions` | A BED file with regions of interest used in QC. By default this is the same as `target_regions`. | `string` |  |  
|  |                                                                                                                               
| `cnv_expected_xy_cn` | A BED file containing expected copy number regions for XY samples. | `string` |  |  |  |                  
| `cnv_expected_xx_cn` | A BED file containing expected copy number regions for XX samples. | `string` |  |  |  |                  
| `cnv_excluded_regions` | A BED file specifying regions to exclude with HiFiCNV or Sawfish, such as centromeres. | `string` |  |  
|  |                                                                                                                               
| `genmod_reduced_penetrance` | A file with gene ids that have reduced penetrance. For use with genmod. | `string` |  |  |  |      
| `genmod_score_config_snvs` | A SNV rank model config file for genmod. | `string` |  |  |  |                                      
| `genmod_score_config_svs` | A SV rank model config file for genmod. | `string` |  |  |  |                                        
| `somalier_sites` | A VCF of known polymorphic sites for somalier | `string` |  |  |  |                                           
| `gens_baf_positions` | Tab-delimited file with variant positions used to calculate B-allele frequencies for Gens inputs. |       
`string` |  |  |  |                                                                                                                
| `gens_panel_of_normals` | Panel of normals (HDF5) used to standardize coverage for Gens inputs                                   
(https://gatk.broadinstitute.org/hc/en-us/articles/360040510031-CreateReadCountPanelOfNormals). This can be generated using        
https://github.com/nf-core/createpanelrefs. For long reads we have used mosdepth to calculate coverage per-bin rather than counting
reads per bin, and then created a panel of normals using the GATK command. | `string` |  |  |  |                                   
| `gens_coverage_bins` | Bed file with bins for which to calculate coverage. Typically 100bp bins across the genome is used in     
Gens. It should be the same as the one used to prepare the panel of normals. | `string` |  |  |  |                                 
| `strdrop_training_set_json` | A JSON file containing the training set for strdrop | `string` |  |  |  |                          
| `peddy_sites` | A file path to a VCF of known polymorphic sites for peddy. You may need to create a custom sites file if you have
incomplete or targeted data. | `string` |  |  |  |                                                                                 
| `sambamba_regions` | A BED file with regions of interest used in sambamba depth. By default this is the same as `qc_regions`. |  
`string` |  |  |  |                                                                                                                
| `alignment_output_format` | Output format for alignment files. Either `bam` or `cram` (accepted: `bam`\|`cram`) | `string` | bam 
|  |  |                                                                                                                            
| `modules_testdata_base_path` | Base URL or local path to location of modules test dataset files | `string` |                     
https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/ |  | True |                                                  
                                                                                                                                   
## Reference genome options                                                                                                        
                                                                                                                                   
Reference genome related files and options required for the workflow.                                                              
                                                                                                                                   
| Parameter | Description | Type | Default | Required | Hidden |                                                                   
|-----------|-----------|-----------|-----------|-----------|-----------|                                                          
| `fasta` | Reference genome | `string` |  |  |  |                                                                                 
                                                                                                                                   
## Institutional config options                                                                                                    
                                                                                                                                   
Parameters used to describe centralised config profiles. These should not be edited.                                               
                                                                                                                                   
| Parameter | Description | Type | Default | Required | Hidden |                                                                   
|-----------|-----------|-----------|-----------|-----------|-----------|                                                          
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |                               
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running        
offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is 
not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this     
parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |                      
| `config_profile_name` | Institutional config name. | `string` |  |  | True |                                                     
| `config_profile_description` | Institutional config description. | `string` |  |  | True |                                       
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |                                   
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |                                                  
| `publish_unannotated_family_svs` | Publish unannotated SVs and CNVs per family and caller. | `boolean` |  |  | True |            
                                                                                                                                   
## Generic options                                                                                                                 
                                                                                                                                   
Less common options for the pipeline, typically set in a config file.                                                              
                                                                                                                                   
| Parameter | Description | Type | Default | Required | Hidden |                                                                   
|-----------|-----------|-----------|-----------|-----------|-----------|                                                          
| `version` | Display version and exit. | `boolean` |  |  | True |                                                                 
| `publish_dir_mode` | Method used to save pipeline results to output directory. (accepted:                                        
`symlink`\|`rellink`\|`link`\|`copy`\|`copyNoFollow`\|`move`) <details><summary>Help</summary><small>The Nextflow `publishDir`     
option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method      
should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for           
details.</small></details>| `string` | copy |  | True |                                                                            
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email
address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit                        
successfully.</small></details>| `string` |  |  | True |                                                                           
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |                                            
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |      
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |                                                  
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging      
service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |                                   
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |                                            
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | 
True |                                                                                                                             
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |  
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |          
| `pipelines_testdata_base_path` | Base URL or local path to location of pipeline test dataset files | `string` |                  
https://raw.githubusercontent.com/genomic-medicine-sweden/test-datasets/6d0a81260c1e8d560eb48f405f405aa0d54c3488/ |  | True |      
| `trace_report_suffix` | Suffix to add to the trace report filename. Default is the date and time in the format                   
yyyy-MM-dd_HH-mm-ss. | `string` |  |  | True |                                                                                     
| `help` | Display the help message. | `['boolean', 'string']` |  |  |  |                                                          
| `help_full` | Display the full detailed help message. | `boolean` |  |  |  |                                                     
| `show_hidden` | Display hidden parameters in the help message (only works when --help or --help_full are provided). | `boolean` |
|  |  |                                                                                                                            
                                                                                                                                   
## Workflow options                                                                                                                
                                                                                                                                   
Workflow options specific to genomic-medicine-sweden/nallo                                                                         
                                                                                                                                   
| Parameter | Description | Type | Default | Required | Hidden |                                                                   
|-----------|-----------|-----------|-----------|-----------|-----------|                                                          
| `bcftools_roh_af_tag` | AF-tag that variants have been annotated with for use in the chromograph subworkflow. By default set to  
`--chromograph_af_tag`. | `string` |  |  |  |                                                                                      
| `chromograph_af_tag` | AF-tag to use in the chromograph subworkflow to generate plots of autozygosity. SNVs need to be annotated 
with allele frequencies specified by this tag for autozygosity plotting to work. | `string` |  |  |  |                             
| `preset` | Enable or disable certain parts of the pipeline by default, depending on data type (`revio`, `pacbio`, `ONT_R10`)     
(accepted: `revio`\|`pacbio`\|`ONT_R10`) | `string` | revio | True |  |                                                            
| `sv_callers` | Which SV callers to use. Several callers can be specified, separated by commas (e.g.                              
sniffles,severus,hificnv,sawfish). The order of the SV callers in this list will determine the priority of the calls when merging  
them if not overwritten by `sv_caller_priority`. | `string` | sniffles,hificnv |  |  |                                             
| `sv_callers_merge_priority` | The order of the SV callers in this list will determine the priority of the calls when merging     
them. All callers that has been specified in `sv_callers` should also be specified here separated by commas (e.g.                  
sniffles,severus,hificnv,sawfish). By default same as `--sv_callers`. | `string` | sniffles,hificnv |  |  |                        
| `sv_callers_to_run` | Which SV callers to run, separated by commas (e.g. sniffles,severus,hificnv,sawfish). By default same as   
`--sv_callers` | `string` | sniffles,hificnv |  |  |                                                                               
| `sv_callers_to_merge` | Which SV callers to merge into family VCFs (that are also used for annotation and ranking), separated by 
commas (e.g. sniffles,severus,hificnv,sawfish). By default same as `--sv_callers` | `string` | sniffles,hificnv |  |  |            
| `snv_caller` | Which short variant software to use (`deepvariant`) (accepted: `deepvariant`) | `string` | deepvariant |  |  |    
| `str_caller` | Which caller to use for short tandem repeat expansions (TRGT or STRdust). (accepted: `trgt`\|`strdust`) | `string`
| trgt |  |  |                                                                                                                     
| `phaser` | Which phasing software to use (`longphase`, `whatshap`, `hiphase`) (accepted: `longphase`\|`whatshap`\|`hiphase`) |   
`string` | longphase |  |  |                                                                                                       
| `hifiasm_mode` | Run hifiasm in hifi-only or hifi-trio mode (`hifi-only`, `trio-binning`) (accepted: `hifi-only`\|`trio-binning`)
| `string` | trio-binning |  |  |                                                                                                  
| `hifiasm_preset` | Hifiasm preset, is set to `--ont` when `--preset ONT_R10` is active. (accepted: ``\|`--ont`) | `string` | None
|  |  |                                                                                                                            
| `run_methbat` | Run methbat for methylation analysis, set to `true` by default when `--preset revio` is active. | `boolean` |  | 
|  |                                                                                                                               
| `run_modkit` | Run modkit for methylation analysis, set to `true` by default when `--preset ONT_R10` is active. | `boolean` |  | 
|  |                                                                                                                               
| `methbat_male_label` | Label used for male samples in methbat profile. | `string` | MALE |  |  |                                 
| `methbat_female_label` | Label used for female samples in methbat profile. | `string` | FEMALE |  |  |                           
| `alignment_processes` | If alignment_processes is bigger than 1, input files will be split and aligned in parallel to reduce     
processing time. | `integer` | 8 |  |  |                                                                                           
| `snv_calling_processes` | If snv_calling_processes is bigger than 1, short variant calling will be done in parallel to reduce    
processing time. | `integer` | 13 |  |  |                                                                                          
| `vep_cache_version` | VEP cache version | `integer` | 110 |  |  |                                                                
| `vep_plugin_files` | Path to a CSV/TSV/JSON/YAML file with vep_files as header, and then paths to vep plugin files. Paths to     
pLI_values.txt and LoFtool_scores.txt are required. | `string` |  |  |  |                                                          
| `force_sawfish_joint_call_single_samples` | Force sawfish to run joint-call on single samples instead of all samples from the    
same family. This effectively causes SVDB to merge the samples into family VCFs instead. | `boolean` |  |  |  |                    
| `filter_variants_hgnc_ids` | A tsv/csv file with a `hgnc_ids` column header, and then one numerical HGNC ID per row. E.g. `4281` 
or `HGNC:4281`. | `string` |  |  |  |                                                                                              
| `filter_snvs_expression` | An expression that is passed to bcftools view to filter SNVs, e.g. --filter_snvs_expression "-e       
'INFO/AQ>60'" | `string` | None |  |  |                                                                                            
| `filter_svs_expression` | An expression that is passed to bcftools view to filter SVs, e.g. --filter_svs_expression "-e          
'INFO/AQ>60'" | `string` | None |  |  |                                                                                            
| `deepvariant_model_type` | Sets the model type used for DeepVariant. This is set automatically using `--preset` by default.      
(accepted: `PACBIO`\|`ONT_R104`) | `string` | PACBIO |  | True |                                                                   
| `minimap2_read_mapping_preset` | Sets the minimap2-preset (-x) for read alignment. This is set automatically using the pipeline  
`--preset` by default. (accepted: `map-hifi`\|`map-ont`\|`lr:hq`\|`lr:hqae`) | `string` | map-hifi |  | True |                     
| `extra_methbat_profile_options` | Extra options to methbat profile. | `string` |  |  | True |                                    
| `extra_modkit_options` | Extra options to modkit, used for test profile. | `string` |  |  | True |                               
| `extra_vep_options` | Extra options to VEP, used for test profile. | `string` |  |  | True |                                     
| `extra_paraphase_options` | Extra options to Paraphase, used for test profile. | `string` |  |  | True |                         
| `extra_hifiasm_options` | Extra options to hifiasm, used for test profile. | `string` |  |  | True |                             
| `extra_yak_options` | Extra options to yak, used for test profile. | `string` |  |  | True |                                     
| `genmod_compound_snv_penalty` | Genmod compound penalty for SNVs. | `integer` | 6 |  |  |                                        
| `genmod_compound_snv_threshold` | Genmod compound threshold for SNVs. | `integer` | 9 |  |  |                                    
| `genmod_compound_sv_penalty` | Genmod compound penalty for SVs. | `integer` | 6 |  |  |                                          
| `genmod_compound_sv_threshold` | Genmod compound threshold for SVs. | `integer` | 9 |  |  |                                      
| `genmod_compound_singleton_snv_penalty` | Genmod compound penalty for SNVs. Applied in families with an affected individual      
without two parents. By default same as `--genmod_compound_snv_penalty`. | `integer` |  |  |  |                                    
| `genmod_compound_singleton_snv_threshold` | Genmod compound threshold for SNVs. Applied in families with an affected individual  
without two parents. By default same as `--genmod_compound_snv_threshold`. | `integer` |  |  |  |                                  
| `genmod_compound_singleton_sv_penalty` | Genmod compound penalty for SVs. Applied in families with an affected individual without
two parents. By default same as `--genmod_compound_sv_penalty`. | `integer` |  |  |  |                                             
| `genmod_compound_singleton_sv_threshold` | Genmod compound threshold for SVs. Applied in families with an affected individual    
without two parents. By default same as `--genmod_compound_sv_threshold`. | `integer` |  |  |  |                                   
| `genmod_compound_trio_snv_penalty` | Genmod compound penalty for SNVs. Applied in families with at least one affected individual 
with two parents. By default same as `--genmod_compound_snv_penalty`. | `integer` |  |  |  |                                       
| `genmod_compound_trio_snv_threshold` | Genmod compound threshold for SNVs. Applied in families with at least one affected        
individual with two parents. By default same as `--genmod_compound_snv_threshold`. | `integer` |  |  |  |                          
| `genmod_compound_trio_sv_penalty` | Genmod compound penalty for SVs. Applied in families with at least one affected individual   
with two parents. By default same as `--genmod_compound_sv_penalty`. | `integer` |  |  |  |                                        
| `genmod_compound_trio_sv_threshold` | Genmod compound threshold for SVs. Applied in families with at least one affected          
individual with two parents. By default same as `--genmod_compound_sv_threshold`. | `integer` |  |  |  |                           
| `plot_chromograph_autozygosity` | Whether to plot chromograph autozygosity plots in the chromograph subworkflow. This requires   
SNVs to be annotated with allele frequencies, and is therefore false by default unless `--bcftools_roh_af_tag` and                 
`--rhocallviz_af_tag` have been set, in which case it is assumed that the user has annotated the SNVs with the correct tag. |      
`boolean` |  |  |  |                                                                                                               
| `plot_chromograph_coverage` | Whether to plot chromograph coverage plots in the chromograph subworkflow. | `boolean` | True |  | 
|                                                                                                                                  
| `pre_vep_snv_filter_expression` | An expression that is passed to bcftools view to filter SNVs before being annotated with VEP,  
e.g. --pre_vep_snv_filter_expression "-e 'INFO/AQ>60'". The expression applies to both the clinical and the research VCFs. |       
`string` | None |  |  |                                                                                                            
| `rhocallviz_af_tag` | AF-tag that variants have been annotated with, for use in the chromograph subworkflow. By default set to   
`--chromograph_af_tag. | `string` |  |  |  |                                                                                       
| `rhocallviz_min_af` | Minimum allele frequency for variants to be included in rhocall viz within the chromograph subworkflow. |  
`number` | 0.0001 |  |  |                                                                                                          
| `rhocallviz_min_qual` | Minimum quality for variants to be included in rhocall viz within the chromograph subworkflow. | `number`
| 10.0 |  |  |                                                                                                                     
| `tiddit_bin_size` | Bin size to use for TIDDIT coverage wig generation in the chromograph subworkflow. | `integer` | 500 |  |  | 
                                                                                                                                   


