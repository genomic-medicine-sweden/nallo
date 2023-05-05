process FIX_TRIO {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    //tuple val(id), val(family_id), val(paternal_id), val(maternal_id), val(sex), val(phenotype), file(bam) 
    tuple val(metameta), path(bam) 

    output:
    path "versions.yml", emit: versions
    tuple stdout, val(metameta), path("*.trio.*"), emit: trio
    //stdout emit: trio
    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in fellen31/skierfe/bin/
    """
    
    PATERNAL_ID=\$(echo $metameta | sed 's/.*paternal_id://' | sed 's/,.*//')
    MATERNAL_ID=\$(echo $metameta | sed 's/.*maternal_id://' | sed 's/,.*//')
    INDIVIDUAL_ID=\$(echo $metameta | sed 's/\\[id://g' | sed 's/,.*//')
    
    echo \$INDIVIDUAL_ID > ID_FILE
    ls *.bam >> BAMS

    IS_EMPTY=0

    if [ -f \$INDIVIDUAL_ID.bam ]; then 
      mv \$INDIVIDUAL_ID.bam \$INDIVIDUAL_ID.kid.trio.bam
      mv \$INDIVIDUAL_ID.bam.bai \$INDIVIDUAL_ID.kid.trio.bam.bai
    else 
      IS_EMPTY=1
      touch x.trio.bam
    fi
    
    if [ -f \$PATERNAL_ID.bam ]; then 
      mv \$PATERNAL_ID.bam \$PATERNAL_ID.dad.trio.bam
      mv \$PATERNAL_ID.bam.bai \$PATERNAL_ID.dad.trio.bam.bai
    else 
      IS_EMPTY=1
      touch x.trio.bam
    fi
    
    
    if [ -f \$MATERNAL_ID.bam ]; then 
      mv \$MATERNAL_ID.bam \$MATERNAL_ID.mom.trio.bam
      mv \$MATERNAL_ID.bam.bai \$MATERNAL_ID.mom.trio.bam.bai
    else 
      IS_EMPTY=1
      touch x.trio.bam
    fi
    
    echo \$IS_EMPTY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
