
process MODKIT_SAMPLEPROBS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "docker.io/fellen31/modkit:v0.5.1-rc1"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), env(probs), emit: probs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"

    """
    #Dynamic determination of 'filter-thresholds'
    modkit \\
        sample-probs \\
        --threads ${task.cpus} \\
        -p 0.1 --interval-size 5000000 --only-mapped \\
        -o . --hist --prefix ${prefix}_modkit_sample-probs \\
        $bam

    probs=\$( awk 'NR>1 {ORS=" "; print "--filter-threshold "\$1":"\$3}' ${prefix}_modkit_sample-probs_thresholds.tsv )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$( modkit --version | sed 's/mod_kit //' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.bedgraph
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$( modkit --version | sed 's/mod_kit //' )
    END_VERSIONS
    """
}
