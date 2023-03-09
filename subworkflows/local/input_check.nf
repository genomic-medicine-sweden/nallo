/*
 * Does NOT check samplesheet, just gets channels from samplesheet
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
  take:
  samplesheet

  main:

  SAMPLESHEET_CHECK ( samplesheet )
    .csv
    .splitCsv ( header:false, sep:',' )
    .map { it -> [ it[0], it[1], it[2]] }
    .set { ch_sample }
  
  emit:
  ch_sample
}
