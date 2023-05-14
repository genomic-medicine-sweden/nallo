process SNIFFLES_MULTISAMPLE {
    tag "sniffles_multisample"
    label 'process_medium'

    conda "bioconda::sniffles=2.0.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:2.0.7--pyhdfd78af_0':
        'quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0' }"

    input:
    path(snfs)

    output:
    path("multisample.sniffles.vcf"), emit: multisample_vcf
    path "versions.yml"                    , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sniffles \\
        --input ${snfs} \\
        --vcf multisample.sniffles.vcf \\
        -t ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --help 2>&1 | grep Version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}
