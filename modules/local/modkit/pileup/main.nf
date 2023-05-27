process MODKIT_PILEUP {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/fellen31/modkit-tabix:latest"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bed.gz"), path("*.bed.gz.tbi"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bgzip <(modkit pileup \\
        $args \\
        --threads ${task.cpus} \\
        --ref $fasta \\
        --log-filepath ${bam.baseName}.modkit.log \\
        ${bam} \\
        -) \\
        -@ {$task.cpus} \\
        -c > ${meta.id}.${bam.baseName}.modkit.bed.gz 

    tabix \\
        -p bed \\
        ${meta.id}.${bam.baseName}.modkit.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$( modkit --version | sed 's/mod_kit //' )
        bgzip: \$( bgzip --version | head -n 1 | sed 's/bgzip (htslib) //g')
        tabix: \$( tabix --version | head -n 1 | sed 's/tabix (htslib) //g')
    END_VERSIONS
    """
}
