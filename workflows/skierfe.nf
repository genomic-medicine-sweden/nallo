////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check mandatory parameters (missing protocol or profile will exit the run.)
if (params.input) { 
    ch_input = file(params.input) 
} else {
    exit 1, 'Input samplesheet not specified!'
}

// Temporary input for extra SNF-files for Sniffles
params.extra_snfs = 'some_value'
params.extra_gvcfs = 'some_value'
if (params.extra_snfs) {
    ch_input_snfs = file(params.extra_snfs)
} else {
    ch_input_snfs = Channel.empty()
  }

// Temporary input for extra gVCF files for GLNexus 

if (params.extra_gvcfs) {
    ch_input_gvcfs = file(params.extra_gvcfs)
} else {
    ch_input_gvcfs = Channel.empty()
  }
////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

/*
 * SUBWORKFLOW:
 */

  include { INPUT_CHECK                      } from '../subworkflows/local/input_check.nf'
  include { INPUT_CHECK as INPUT_SNFS_CHECK  } from '../subworkflows/local/input_check.nf'
  include { INPUT_CHECK as INPUT_GVCFS_CHECK } from '../subworkflows/local/input_check.nf'
  include { STRUCTURAL_VARIANT_CALLING       } from '../subworkflows/local/structural_variant_calling'
  include { SNP_CALLING                      } from '../subworkflows/local/snp_calling'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow SKIERFE {

  /*
   * SUBWORKFLOW: Read in samplesheet, validate and stage input files
  */

  // Need to set a channel with sample names
  
  INPUT_CHECK ( ch_input )
    .set { ch_sample }


// TODO: Fix this very ugly fix
if (params.extra_snfs != 'some_value') {
  INPUT_SNFS_CHECK ( ch_input_snfs )
    .set { ch_extra_snfs }
} else {
    ch_extra_snfs = Channel.empty()
  }

if (params.extra_snfs != 'some_value') {
  INPUT_GVCFS_CHECK ( ch_input_gvcfs )
    .set { ch_extra_gvcfs }
} else {
    ch_extra_gvcfs = Channel.empty()
  }
  /*
   * SUBWORKFLOW: Structural variant calling
   */

  // Take this from sampleheet instead of, but could maybe check automatically for index?:
  //bam = channel.fromFilePairs("$baseDir/data/raw/*.{bam,bai}", flat:true, checkIfExists:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
  
  ch_fasta = channel.fromPath(params.reference_fasta)
  ch_fai = channel.fromPath(params.reference_index)

  STRUCTURAL_VARIANT_CALLING ( ch_sample, ch_extra_snfs )
  SNP_CALLING ( ch_sample, ch_fasta, ch_fai, ch_extra_gvcfs )
  
}
