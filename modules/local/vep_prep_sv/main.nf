process VEP_PREP_SV {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/99ae4e6605f4b9f9d3197dcc856f867fbbd1e6a334dee63957e56f908ed4599d/data'
        : 'community.wave.seqera.io/library/htslib_python:ef76279c6a229e2e'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('vep_prep_sv'), eval("vep_prep_sv.py --version"), emit: versions_vep_prep_sv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    vep_prep_sv.py ${args} ${vcf} > ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
