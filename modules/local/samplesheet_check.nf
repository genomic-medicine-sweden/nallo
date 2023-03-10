process SAMPLESHEET_CHECK {
  tag "$samplesheet"
  time '15min'
  label 'small'
  input:
  path samplesheet

  output:
  path '*.csv', emit: csv

  script:
  """
  mv ${samplesheet} ${samplesheet.getSimpleName()}.unchecked.csv 
  """
  

  }
