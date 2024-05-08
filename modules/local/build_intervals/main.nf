// https://github.com/nf-core/sarek/blob/59026dc07633edb83aab3bfb2f65f79db38437a1/modules/local/build_intervals/main.nf

process BUILD_INTERVALS {
    tag "$meta.id"

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(fasta_fai)

    output:
    tuple val(meta), path("${fasta_fai.baseName}.bed") , emit: bed
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fasta_fai} > ${fasta_fai.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
