process COLLAPSE_MITOCHONDRIAL_GT {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbeb20c898a76bec809629320ece9e1f84a3e355e96568bfbe14b9f411bdf3e7/data'
        : 'community.wave.seqera.io/library/htslib_python:9c6265e98ef06930'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('collapse_mitochondrial_gt'), eval("collapse_mitochondrial_gt.py --version"), emit: versions_collapse_mitochondrial_gt, topic: versions
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), emit: versions_bgzip, topic: versions
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), emit: versions_tabix, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    collapse_mitochondrial_gt.py ${args} ${vcf} > ${prefix}.vcf
    bgzip --threads ${task.cpus} ${prefix}.vcf
    tabix --threads ${task.cpus} -p vcf ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | bgzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
