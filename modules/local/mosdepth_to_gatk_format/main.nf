process MOSDEPTH_TO_GATK_FORMAT {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"
    
    input:
    tuple val(meta), path(mosdepth_file)
    path(mosdepth_header_template)

    output:
    tuple val(meta), path("*.tsv"), emit: output
    path "versions.yml", emit: versions

    script:
    """
    # Insert sample name into header template
    sed "s/SAMPLE_NAME/${meta.id}/" ${mosdepth_header_template} > ${meta.id}.header

    # Round mosdepth output to nearest integer
    # Position starting at 1
    gzip -cd $mosdepth_file |\
        awk '
            BEGIN {OFS="\t"}
            { 
            \$4 = int(\$4 + 0.5)
            \$2++
            print \$1, \$2, \$3, \$4}
        ' > ${meta.id}.body

    # Skip special chromosomes. Gens only accepts chr1-22, chrX, chrY, chrM
    # grep -vP "^\w+_|^chrE" ${meta.id}_all.body > ${meta.id}_target.body
    # grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\t' ${meta.id}_all.body > ${meta.id}_target.body

    cat ${meta.id}.header ${meta.id}.body > ${meta.id}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | head -n 1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | head -n 1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """
}