process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/google/deepvariant:1.5.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(bam), path(bai), val(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.deepvariant.vcf.gz")  , emit: vcf // Need something to distinguish from gVCF
    tuple val(meta), path("*.g.vcf.gz"), emit: gvcf
    tuple val(meta), path("*.g.vcf.gz"), val(intervals), emit: gvcf_region
    path "*.html"                      , emit: html
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def model_type = task.ext.model_type ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def regions    = intervals ? "--regions ${intervals}" : ""

    def underscore_regions = intervals.replaceAll(/[:-]/, "_")
    def output_name = intervals ? "${prefix}" + "." + "${underscore_regions}" : "${prefix}"

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --ref=${fasta} \\
        --reads=${bam} \\
        --output_vcf=${output_name}.deepvariant.vcf.gz \\
        --output_gvcf=${output_name}.g.vcf.gz \\
        --sample_name=${meta.id} \\
        --num_shards=${task.cpus} \\
        --model_type=$model_type \\
        ${regions} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}

