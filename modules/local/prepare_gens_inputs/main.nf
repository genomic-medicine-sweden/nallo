// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.

process PREPAREGENSINPUTDATA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/902de3e3c7bf48ca1ece3dfc74eec40431a264c9ca50c182910702e1518fe265/data':
        'community.wave.seqera.io/library/pip_gens-input-data-tools:ad7b7f1a90b0d4ca' }"

    input:
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(read_counts_gz), path(gvcf), path(gvcf_tbi)
    path baf_positions

    output:
    tuple val(meta), path("*.cov.bed.gz")     , emit: cov_gz
    tuple val(meta), path("*.cov.bed.gz.tbi") , emit: cov_tbi
    tuple val(meta), path("*.baf.bed.gz")     , emit: baf_gz
    tuple val(meta), path("*.baf.bed.gz.tbi") , emit: baf_tbi
    path "versions.yml"                       , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "The gens pre-processing module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    new_gens_command.py \\
        --coverage $read_counts_gz \\
        --gvcf $gvcf \\
        --label $prefix \\
        --baf_positions $baf_positions \\
        $args \\
        --outdir \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preparegensinputdata: \$(new_gens_command.py --version | sed 's/^v//')
    END_VERSIONS
    """
    // FIXME: Check if the version is clean (i.e. trim v)

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.cov.bed.gz
    touch ${prefix}.cov.bed.gz.tbi
    echo "" | gzip > ${prefix}.baf.bed.gz
    touch ${prefix}.baf.bed.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preparegensinputdata: \$(new_gens_command.py --version | sed 's/^v//')
    END_VERSIONS
    """
}
