#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process sniffles {
  cpus = 36
  container = 'docker.io/library/sniffles:2.0.7'
  publishDir 'data/interim/sv-calling/sniffles/', mode: 'copy'
  
  input:
    tuple val(sampleId), path(reads), path(index)
  output:
    path '*.vcf', emit: samples_vcf
    path '*.snf', emit: samples_snf

  """
  sniffles \
    --input $reads \
    --vcf ${sampleId}.sniffles.ng.vcf \
    --snf ${sampleId}.sniffles.ng.snf \
    --threads ${task.cpus} \
    --non-germline
  """
}

process sniffles_combined_calling {
  cpus = 36
  container = 'docker.io/library/sniffles:2.0.7'
  publishDir 'data/interim/sv-calling/sniffles/', mode: 'copy'
  
  input:
    path(snf)
  output:
    path '*.vcf'

  """
  sniffles --input ${snf} --vcf multisample.ng.vcf --threads ${task.cpus}
  """

  }

workflow {
  bam = channel.fromFilePairs("$baseDir/data/alignments/*.{bam,bai}", flat:true, checkIfExists:true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
  sniffles(bam)
  sniffles_combined_calling(sniffles.out.samples_snf.collect())
}
