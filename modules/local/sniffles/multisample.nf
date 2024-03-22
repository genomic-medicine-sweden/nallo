process SNIFFLES_MULTISAMPLE {
    tag "sniffles_multisample"
    label 'process_medium'

    conda "bioconda::sniffles=2.0.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sniffles:2.0.7--pyhdfd78af_0':
        'quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0' }"

    input:
    path(snfs)
    tuple val(meta), path(reference)
    path(tandem_file)

    output:
    path("multisample.sniffles.vcf"), emit: vcf
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tandem_repeats = tandem_file ? "--tandem-repeats ${tandem_file}" : ''
    """
    sniffles \\
        --input ${snfs} \\
        --vcf multisample.sniffles.vcf \\
        -t ${task.cpus} \\
        --reference $reference \\
        $tandem_repeats \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles_multisample: \$(sniffles --help 2>&1 | grep Version |sed 's/^.*Version: //')
    END_VERSIONS
    """
}
