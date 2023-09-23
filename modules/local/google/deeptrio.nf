process DEEPTRIO {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/google/deepvariant:deeptrio-1.4.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(individual_bam), path(individual_bai), val(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta), path(paternal_bam), path(paternal_bai)
    tuple val(meta), path(maternal_bam), path(maternal_bai)

    output:
    tuple val(meta), path("*.{kid,dad,mom}.vcf.gz"),   emit: vcf
    tuple val(meta), path("*.{kid,dad,mom}.g.vcf.gz"), emit: gvcf
    path "versions.yml"                        ,  emit: versions
    path "*.html", emit: html

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def model_type = task.ext.model_type ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    //def regions = intervals ? "--regions ${intervals}" : ""

    """
    /opt/deepvariant/bin/deeptrio/run_deeptrio \\
        --ref=${fasta} \\
        --reads_child=${individual_bam} \\
        --reads_parent1=${paternal_bam} \\
        --reads_parent2=${maternal_bam} \\
        --output_vcf_child=${meta.id}.${meta.id}.kid.vcf.gz \\
        --output_vcf_parent1=${meta.id}.${meta.paternal_id}.dad.vcf.gz \\
        --output_vcf_parent2=${meta.id}.${meta.maternal_id}.mom.vcf.gz \\
        --output_gvcf_child=${meta.id}.${meta.id}.kid.g.vcf.gz \\
        --output_gvcf_parent1=${meta.id}.${meta.paternal_id}.dad.g.vcf.gz \\
        --output_gvcf_parent2=${meta.id}.${meta.maternal_id}.mom.g.vcf.gz \\
        --sample_name_child="${meta.id}.${meta.id}" \\
        --sample_name_parent1="${meta.id}.${meta.paternal_id}" \\
        --sample_name_parent2="${meta.id}.${meta.maternal_id}" \\
        --num_shards=${task.cpus} \\
        --model_type=$model_type \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptrio: \$(echo \$(/opt/deepvariant/bin/deeptrio/run_deeptrio --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}

