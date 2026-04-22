process GATK4_CLEANSAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(bam)
    val(create_index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val("${task.process}"), val('gatk'), eval("gatk CleanSam --version | grep GATK | sed 's/.*(GATK) v//'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def index = create_index? "--CREATE_INDEX true" : ""

    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    gatk \\
        CleanSam \\
        ${args} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.bam \\
        ${index} \\
        --use_jdk_inflater true \\
        --use_jdk_deflater true
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index = create_index? "touch ${prefix}.bam.bai" : ""
    if ("${bam}" == "${prefix}.bam") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.bam
    ${index}
    """
}
