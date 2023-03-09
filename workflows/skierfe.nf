////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check mandatory parameters (missing protocol or profile will exit the run.)
if (params.input) { 
    ch_input = file(params.input) 
} else {
    exit 1, 'Input samplesheet not specified!'
}

params.input_snfs = 'some_value'

if (params.input_snfs) {
    ch_input_snfs = file(params.extra_snfs)
} else {
    ch_input_snfs = Channel.empty()
  }

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

/*
 * SUBWORKFLOW:
 */

  include { INPUT_CHECK                      } from '../subworkflows/local/input_check.nf'
  include { INPUT_CHECK as INPUT_SNFS_CHECK  } from '../subworkflows/local/input_check.nf'
  include { STRUCTURAL_VARIANT_CALLING       } from '../subworkflows/local/structural_variant_calling'

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

  INPUT_SNFS_CHECK ( ch_input_snfs )
    .set { ch_extra_snfs }
  
  /*
   * SUBWORKFLOW: Structural variant calling
   */

  // How to pass along metadata?
  
  // Take this from sampleheet instead of:
  //bam = channel.fromFilePairs("$baseDir/data/raw/*.{bam,bai}", flat:true, checkIfExists:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
  
  STRUCTURAL_VARIANT_CALLING ( ch_sample, ch_extra_snfs )
  

}
