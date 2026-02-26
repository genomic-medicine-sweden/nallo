process PORTELLO {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/fellen31/portello:0.7.02"
    input:
    tuple val(meta), path(asm_to_ref_bam), path(asm_to_ref_bai), path(read_to_asm_bam), path(read_to_asm_bai)
    tuple val(meta2), path(fasta)
    val bam_index_extension

    output:
    tuple val(meta), path("*_remapped.bam"), emit: bam
    tuple val(meta), path("*_remapped.bam.${bam_index_extension}"), emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    """
    portello \
        --threads ${task.cpus} \
        --ref $fasta \
        --assembly-to-ref $asm_to_ref_bam \
        --read-to-assembly $read_to_asm_bam \
        --input-assembly-mode partially-phased \
        --remapped-read-output - \
        --unassembled-read-output unassembled.bam |\
        samtools sort -@${task.cpus} - --write-index -o ${prefix}_remapped_no_rg.bam##idx##${prefix}_remapped_no_rg.bam.${bam_index_extension}

    samtools addreplacerg -@ ${task.cpus} -o ${prefix}_remapped.bam ${prefix}_remapped_no_rg.bam -r @RG\\\\tID:${meta.id}\\\\tSM:${meta.id}
    samtools index -@ ${task.cpus} ${prefix}_remapped.bam
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch remapped.sort.bam
    """
}
