process GLNEXUS {
    tag "glnexus_mutlisample"
    label 'process_high'

    container 'quay.io/mlin/glnexus:v1.2.7'
    
     // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DEEPVARIANT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    
    publishDir 'data/interim/short_variant_calling/glnexus'
    
    input:
    path(gvcfs)

    output:
    path("*.vcf.gz"),  emit: multisample_bcf
    path "versions.yml",  emit: versions
    
    """
    glnexus_cli --config DeepVariant_unfiltered --threads ${task.cpus} ${gvcfs} | bcftools view - \\
      | bgzip -@ ${task.cpus} -c \\
      > multisample.vcf.gz
     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: \$(echo \$(glnexus_cli 2>&1) |  head -n 1 | sed 's/^.*glnexus_cli release v//; s/ .*\$//' )
    END_VERSIONS
    """
}
