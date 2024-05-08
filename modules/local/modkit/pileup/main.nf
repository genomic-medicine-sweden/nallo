process MODKIT_PILEUP {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/fellen31/modkit:0.2.5"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(bed)

    output:
    tuple val(meta), path("${meta.id}.bed"), emit: bed, optional: true
    // This will break if there are other tags than 1 and 2 I guess
    tuple val(meta), path("${meta.id}.1.bed")         , emit: haplotype_1, optional: true
    tuple val(meta), path("${meta.id}.2.bed")         , emit: haplotype_2, optional: true
    tuple val(meta), path("${meta.id}.ungrouped.bed") , emit: ungrouped  , optional: true
    path "*.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def include_bed = bed ? "--include-bed ${bed}" : ''
    """
    modkit pileup \\
        $args \\
        $include_bed \\
        --threads ${task.cpus} \\
        --ref $fasta \\
        --log-filepath ${bam.baseName}.modkit.log \\
        ${bam} \\
        ${meta.id}.${bam.baseName}.bed

    if test -d ${meta.id}.${bam.baseName}.bed; then
        if test -f ${meta.id}.${bam.baseName}.bed/1.bed; then mv ${meta.id}.${bam.baseName}.bed/1.bed ${meta.id}.1.bed; fi
        if test -f ${meta.id}.${bam.baseName}.bed/2.bed; then mv ${meta.id}.${bam.baseName}.bed/2.bed ${meta.id}.2.bed; fi
        if test -f ${meta.id}.${bam.baseName}.bed/ungrouped.bed; then mv ${meta.id}.${bam.baseName}.bed/ungrouped.bed ${meta.id}.ungrouped.bed; fi
    else
        mv ${meta.id}.${bam.baseName}.bed ${meta.id}.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$( modkit --version | sed 's/mod_kit //' )
    END_VERSIONS
    """
}
