process COLLAPSE_MITOCHONDRIAL_GT {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/99ae4e6605f4b9f9d3197dcc856f867fbbd1e6a334dee63957e56f908ed4599d/data'
        : 'community.wave.seqera.io/library/htslib_python:ef76279c6a229e2e'}"

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
