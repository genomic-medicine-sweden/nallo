process SNIFFLES_SINGLE_SAMPLE {

  tag "$sampleId"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
  'https://depot.galaxyproject.org/singularity/sniffles:2.0.7--pyhdfd78af_0' :
  'quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0' }"
  
  publishDir 'data/interim/sv-calling/sniffles/', mode: 'copy'
  cpus = 8
  time '1h'

  input:
    tuple val(sampleId), path(reads), path(index)
  
  output:
    tuple val(sampleId), path("*.sniffles.vcf"), emit: sv_vcf
    tuple val(sampleId), path("*.sniffles.snf"), emit: sv_snf

  """
  sniffles \
    --input $reads \
    --vcf ${sampleId}.sniffles.vcf \
    --snf ${sampleId}.sniffles.snf \
    --threads ${task.cpus} 
  """
}

