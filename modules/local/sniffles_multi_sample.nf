
process SNIFFLES_MULTI_SAMPLE {
  
  tag "$sampleId"
  container 'quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0'
  publishDir 'data/interim/sv-calling/sniffles/', mode: 'copy'
  cpus = 16

  input:
      path(snfs)
//    tuple val(sampleId), path(extra_snfs)
  output:
    path("*.sniffles.vcf"), emit: multisample_vcf

  """
  sniffles \
    --input ${snfs} \
    --vcf multisample.sniffles.vcf \
    --threads ${task.cpus}
  """

}
