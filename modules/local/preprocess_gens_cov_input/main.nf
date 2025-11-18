process PREPROCESS_GENS_COV_INPUT {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"
    
    input:
    tuple val(meta), path(mosdepth_file)
    path(mosdepth_header_template)

    output:
    tuple val(meta), path("${meta.id}.tsv"), emit: output
    path "versions.yml", emit: versions

    script:
    // FIXME: Should this be done in a Python util?
    // zcat $mosdepth_file | \
    //     awk '{ 
    //         val = $4 +0
    //         rounded = (val - int(val) >= 0.5) ? int(val) + 1 :int(val)
    //         $4 = rounded
    //         $2 = $2 + 1
    //         OFS = "\t"
    //         print $1, $2, $3, $4
    //     }' > ${meta.id}.tmp
    """
    
    echo "Hi" > ${meta.id}.tmp

    sed "s/SAMPLE_NAME/${meta.id}/" ${mosdepth_header_template} > ${meta.id}.sample
    cat ${meta.id}.sample > ${meta.id}.tsv
    cat ${meta.id}.tmp >> ${meta.id}.tsv

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