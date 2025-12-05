process PBCPGTOOLS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pb-cpg-tools:3.0.0--h9ee0642_0"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)
    val(tool)

    output:
    tuple val(meta), path("*.combined.bed.gz"), path("*.combined.bw"),  emit: bed_and_bw
    //tuple val(meta), path("*.combined.bw"),     emit: bigwig
    tuple val(meta), path("*.tsv"),             emit: tsv,      optional: true
    tuple val(meta), path("*.txt"),             emit: txt,      optional: true
    path "versions.yml",                        emit: versions

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
                --output-prefix ${prefix}
            ;;
        "cpg_pileup")
            cpg_pileup \\
                --bam ${bam} \\
                --ref ${reference} \\
                --output-prefix ${prefix}
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
    touch ${prefix}.bed.gz
    touch ${prefix}.bw
    touch ${prefix}.tsv
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pb-cpg-tools: 2.3.2
    END_VERSIONS
    """
}