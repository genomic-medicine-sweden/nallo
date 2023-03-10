
process SNIFFLES_MULTI_SAMPLE {
  
  tag "$sampleId"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
  'https://depot.galaxyproject.org/singularity/sniffles:2.0.7--pyhdfd78af_0' :
  'quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0' }"

  publishDir 'data/interim/sv-calling/sniffles/', mode: 'copy'
  cpus = 16
  time '24h'

  input:
  path(snfs)
  output:
  path("*.sniffles.vcf"), emit: multisample_vcf

  """
  sniffles \
    --input ${snfs} \
    --vcf multisample.sniffles.vcf \
    --threads ${task.cpus}
  """

}
