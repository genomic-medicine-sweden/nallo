process GENERATE_MOSDEPTH_GATK_HEADER {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"
    
    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.mosdepth_gatk_header.tsv"), emit: header
    path "versions.yml", emit: versions

    script:
    """
    samtools view -H ${bam} > bam_header.txt
    printf "@RG\\tID:GATKCopyNumber\\tSM:${meta.id}\\n" >> bam_header.txt
    printf "CONTIG\\tSTART\\tEND\\tCOUNT\\n" >> bam_header.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.mosdepth_gatk_header.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}