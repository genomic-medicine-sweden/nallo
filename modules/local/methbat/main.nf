process METHBAT_PROFILE {

    tag "$prefix"
    label 'process_medium'    

    conda "bioconda::methbat=0.16.1"
    container "quay.io/biocontainers-methbat:0.17.0--h9ee0642_0"

    input:
    tuple val(prefix), path(bedgz), path(bw), path(regions)

    output:
    path("${prefix}.meth_region_stats.tsv"), emit: profile
    path "versions.yml", emit: versions

    script:
    def region_cmd = regions ? "--input-regions ${regions}" : ""

    """
    mkdir -p ${prefix}

    cp ${bedgz} ${prefix}/${prefix}.combined.bed.gz
    gunzip -f ${prefix}/${prefix}.combined.bed.gz

    cp ${bw} ${prefix}/${prefix}.combined.bw
    
    ${regions ? "cp ${regions} ${prefix}/regions.bed" : ""}

    # Run MethBat without cd; use full paths
    methbat profile \
        --input-prefix ${prefix}/${prefix} \
        ${regions ? "--input-regions ${prefix}/regions.bed" : ""} \
        --output-region-profile ${prefix}.meth_region_stats.tsv \
        --profile-label ALL
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methbat: \$(methbat --version 2>&1 | grep -o 'MethBat [0-9.]*' | cut -d' ' -f2 || echo "unknown")
    END_VERSIONS
    """
}



process METHBAT_COMPARE {
    tag "cohort_methbat"
    label 'process_medium'

    conda "bioconda::methbat=0.16.1"
    container "quay.io/biocontainers/methbat:0.16.1--h9ee0642_0"

    input:
    path cohort_profiles      // list/array of paths
    val baseline_category
    val compare_category

    output:
    path "comparison_results.tsv", emit: comparison
    path "versions.yml",          emit: versions

    script:

    """
    methbat compare \
        --input-profile ${cohort_profiles.join(' ')} \
        --output-comparison comparison_results.tsv \
        --baseline-category ${baseline_category} \
        --compare-category ${compare_category}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methbat: \$(methbat --version 2>&1 | grep -o 'MethBat [0-9.]*' | cut -d' ' -f2 || echo "0.16.1")
    END_VERSIONS
    """
}
