process GENERATE_GENS_DATA {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/ae6d3694118aeb46be476fd3c59a026a177a588754039effa7ad76c5f47713f6/data':
        'community.wave.seqera.io/library/pbgzip_tabix_python:04040cdab96d0c32' }"

    input:
    tuple val(meta), path(coverage), path(gvcf), path(gvcf_tbi)
    path(baf_positions)

    output:
    tuple val(meta), path("*cov.bed.gz"), path("*cov.bed.gz.tbi"), emit: cov_bed_tbi
    tuple val(meta), path("*baf.bed.gz"), path("*baf.bed.gz.tbi"), emit: baf_bed_tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    generate_gens_data.py \\
        --label ${prefix} \\
        --coverage ${coverage} \\
        --gvcf ${gvcf} \\
        --baf_positions ${baf_positions} \\
        --bgzip_tabix_output \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_gens_data: \$(echo \$(generate_gens_data.py --version) )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.cov.bed.gz
    touch ${prefix}.cov.bed.gz.tbi
    touch ${prefix}.baf.bed.gz
    touch ${prefix}.baf.bed.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_gens_data: \$(echo \$(generate_gens_data.py --version) )
    END_VERSIONS
    """
}
