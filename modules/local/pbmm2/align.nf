// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process PBMM2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::pbmm2=1.10.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmm2:1.10.0--h9ee0642_0':
        'quay.io/biocontainers/pbmm2:1.10.0--h9ee0642_0'}"

    publishDir 'data/interim/aligned_reads/pbmm2'

    input:
    tuple val(meta), path(fastq), path(mmi)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam_bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbmm2 align ${mmi} ${fastq} ${meta.id}.bam --preset CCS --sort -j ${task.cpus} -J 4

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$(pbmm2 --version 2>&1 | grep pbmm2 | head -n 1 | sed 's/^.*pbmm2 //; s/ .*\$//')
    END_VERSIONS
    """
}
