process WHATSHAP_HAPLOTAG {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::whatshap=1.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:1.7--py310h30d9df9_0' :
        'quay.io/biocontainers/whatshap:1.7--py310h30d9df9_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(bam), path(bai)
    tuple val(meta3), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // TODO: Include samtools in container and pipe to samtools view instead
    """
    whatshap haplotag \\
        $args \\
        -o ${bam.baseName}.haplotagged.bam \\
        --reference $fasta \\
        --output-threads $task.cpus \\
        ${vcf} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
    END_VERSIONS
    """
}
