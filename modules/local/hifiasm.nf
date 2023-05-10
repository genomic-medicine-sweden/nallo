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

process HIFIASM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::hifiasm=0.19.5"
    container 'quay.io/biocontainers/hifiasm:0.19.5--h5b5514e_0'

    publishDir 'data/interim/hifiasm/'


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.asm.bp.hap1.p_ctg.gfa"), path("${meta.id}.asm.bp.hap2.p_ctg.gfa"), emit: phased_contigs
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hifiasm \\
        -o ${prefix}.asm \\
        -t ${task.cpus} \\
        ${reads} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiasm: \$(hifiasm -v 2>&1)
    END_VERSIONS
    """
}
