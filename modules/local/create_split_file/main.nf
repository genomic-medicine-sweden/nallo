process CREATE_SPLIT_FILE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), val(sample_ids)
    val suffix

    output:
    tuple val(meta), path("*.txt"), emit: txt

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outfile_text = sample_ids.collect { "$it\\t-\\t$it$suffix" }.join("\\n")
    """
    echo -e "$outfile_text" > ${prefix}.txt
    """
    // Not adding stub because the output is used for joining by Nextflow and would fail if missing
}
