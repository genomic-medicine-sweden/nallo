process PBCPGTOOLS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pb-cpg-tools:2.3.2--pyhdfd78af_0':
        'biocontainers/pb-cpg-tools:2.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)
    val(tool)

    output:
    tuple val(meta), path("*.bed"),     emit: bed,     optional: true
    tuple val(meta), path("*.bigwig"),  emit: bigwig,  optional: true
    tuple val(meta), path("*.tsv"),     emit: tsv,     optional: true
    tuple val(meta), path("*.txt"),     emit: txt,     optional: true
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    case "${tool}" in
        "aligned_bam_to_cpg_scores")
            aligned_bam_to_cpg_scores \\
                --bam ${bam} \\
                --ref ${reference} \\
                --output-prefix ${prefix} \\
                ${args}
            ;;
        "cpg_pileup")
            cpg_pileup \\
                --bam ${bam} \\
                --ref ${reference} \\
                --output-prefix ${prefix} \\
                ${args}
            ;;
        "calculate_pmd")
            calculate_pmd \\
                --bam ${bam} \\
                --ref ${reference} \\
                --output-prefix ${prefix} \\
                ${args}
            ;;
        *)
            echo "Unknown pb-cpg-tools command: ${tool}"
            exit 1
            ;;
    esac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pb-cpg-tools: \$(aligned_bam_to_cpg_scores --version 2>&1 | grep -o 'pb-cpg-tools [0-9.]*' | cut -d' ' -f2 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.bigwig
    touch ${prefix}.tsv
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pb-cpg-tools: 2.3.2
    END_VERSIONS
    """
}