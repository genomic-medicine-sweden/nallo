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
        'oras://community.wave.seqera.io/library/pip_gens-input-data-tools:97f213f547c014d1':
        'community.wave.seqera.io/library/pip_gens-input-data-tools:ad7b7f1a90b0d4ca' }"

    input:
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(gvcf), path(gvcf_tbi), path(read_counts_gz)
    path baf_positions

    output:
    tuple val(meta), path("*.cov.bed.gz")     , emit: cov_gz
    tuple val(meta), path("*.cov.bed.gz.tbi") , emit: cov_tbi
    tuple val(meta), path("*.baf.bed.gz")     , emit: baf_gz
    tuple val(meta), path("*.baf.bed.gz.tbi") , emit: baf_tbi
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "The gens pre-processing module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_gens_data.py \\
        --coverage $read_counts_gz \\
        --gvcf $gvcf \\
        --label $prefix \\
        --baf_positions $baf_positions \\
        $args \\
        --outdir \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preparegensinputdata: \$(generate_gens_data.py --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cov.bed.gz
    touch ${prefix}.cov.bed.gz.tbi
    touch ${prefix}.baf.bed.gz
    touch ${prefix}.baf.bed.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preparegensinputdata: \$(generate_gens_data.py --version)
    END_VERSIONS
    """
}
