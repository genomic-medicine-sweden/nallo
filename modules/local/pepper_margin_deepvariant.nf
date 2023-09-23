process PEPPER_MARGIN_DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/kishwars/pepper_deepvariant:r0.8"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(bam), path(bai), val(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")  , emit: vcf
    tuple val(meta), path("${prefix}.g.vcf.gz"), emit: gvcf
    path "*.html"                              , emit: html
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def model_type = task.ext.model_type ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    //def regions = intervals ? "--regions ${intervals}" : ""

    """
    run_pepper_margin_deepvariant call_variant \\
        -f ${fasta} \\
        -b ${bam} \\
        -o . \\
        -p ${meta.id} \\
        -s ${meta.id} \\
        -t ${task.cpus} \\
        $model_type \\
        --gvcf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}

