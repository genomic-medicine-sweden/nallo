process PORTELLO {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/fellen31/portello:0.7.02"
    input:
    tuple val(meta), path(asm_to_ref_bam), path(asm_to_ref_bai), path(read_to_asm_bam), path(read_to_asm_bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("remapped.sort.bam"), emit: bam
    tuple val(meta), path("remapped.sort.bam.*i"), emit: index

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
        samtools sort -@${task.cpus} - --write-index -o remapped.sort.bam
    samtools index remapped.sort.bam
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch remapped.sort.bam
    """
}
