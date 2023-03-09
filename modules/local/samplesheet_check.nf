process SAMPLESHEET_CHECK {
  tag "$samplesheet"

  input:
  path samplesheet

  output:
  path '*.csv', emit: csv

  script:
  """
  mv ${samplesheet} ${samplesheet.getSimpleName()}.unchecked.csv 
  """
  

  }
