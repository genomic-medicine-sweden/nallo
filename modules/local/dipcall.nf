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

process DIPCALL {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::dipcall=0.3"
    container 'quay.io/biocontainers/dipcall:0.3--0'

    publishDir 'data/interim/dipcall/'


    input:
    tuple val(meta), path(haplotype_1), path(haplotype_2), path(reference), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: sv_calls
    tuple val(meta), path("*.bed"), emit: confident_regions
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    run-dipcall -t ${task.cpus} $args ${prefix} ${reference} ${haplotype_1} ${haplotype_2} > ${prefix}.mak
    
    make -j2 -f ${prefix}.mak
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dipcall: \$(dipcall-aux.js version)
        make: \$(make -v |head -1 | sed 's/.* //')
    END_VERSIONS
    """
}
