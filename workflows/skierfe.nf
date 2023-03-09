////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

/*
 * SUBWORKFLOW:
 */

 include { STRUCTURAL_VARIANT_CALLING       } from '../subworkflows/local/structural_variant_calling'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow SKIERFE {

  /*
   * SUBWORKFLOW: Read in samplesheet, validate and stage input files
  */

  // Need to set a channel with sample names
 
  /*
   * SUBWORKFLOW: Structural variant calling
   */

  // How to pass along metadata?
  
  // Take this from sampleheet instead of:
  // bam = channel.fromFilePairs("$baseDir/data/alignments/*.{bam,bai}", flat:true, checkIfExists:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }

  STRUCTURAL_VARIANT_CALLING ( ch_view_bam )

}
