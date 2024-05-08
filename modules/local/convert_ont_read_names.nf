process CONVERT_ONT_READ_NAMES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input), path(index)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"),  emit: bam_bai
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt bam") ? "bam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    input.getExtension()
    if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    // The SED looks for readnames in a BAM file
    // and replaces all hyphens - with underscores _,
    // then it adds /77923358/ccs/ to the read name
    // to make ONT read names look like they come from PacBio.
    //
    // The first part of the sed command: '/^[^@]/s/-/_/g
    // Matches lines that do not start with @ (matches the reads) and replaces all occurrences of - with _.
    //
    // The second part: ;/^[^@]/s/^([^[:space:]]+)/\\1\\/77923358\\/ccs/'
    // Also matches lines that do not start with @,
    // and replaces the first sequence of non-space characters (readname) with itself (\\1), followed by /77923358/ccs/.
    """
    samtools view -x MM,ML --threads ${(task.cpus-1)/2} -h $input |\\
    sed -E '/^[^@]/s/-/_/g;/^[^@]/s/^([^[:space:]]+)/\\1\\/77923358\\/ccs/' |\\
    samtools \\
        view \\
        --threads ${(task.cpus-1)/2} \\
        $args \\
        -o ${prefix}.${file_type}

    samtools index -@ ${task.cpus-1} ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
