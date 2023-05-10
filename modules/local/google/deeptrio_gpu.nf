process DEEPTRIO_GPU {
    tag "$meta.id"
    label 'process_gpu'
    
    container "google/deepvariant:deeptrio-1.4.0-gpu"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    publishDir 'data/interim/short_variant_calling/deeptrio'

    input:
    tuple val(meta), path(fasta), path(fai), path(files)

    output:
    tuple val(meta), path("*.{kid,dad,mom}.vcf.gz"),  emit: vcf
    tuple val(meta), path("*.{kid,dad,mom}.g.vcf.gz"),  emit: gvcf
    path "*.html", emit: html
    path "versions.yml"                        ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    //def regions = intervals ? "--regions ${intervals}" : ""

    """
    /opt/deepvariant/bin/deeptrio/run_deeptrio \\
        --ref=${fasta} \\
        --reads_child=${meta.id}.kid.trio.bam \\
        --reads_parent1=${meta.paternal_id}.dad.trio.bam \\
        --reads_parent2=${meta.maternal_id}.mom.trio.bam \\
        --model_type=PACBIO \\
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
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptrio: \$(echo \$(/opt/deepvariant/bin/deeptrio/run_deeptrio --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}

