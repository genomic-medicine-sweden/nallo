process CREATE_PEDIGREE_FILE {
    tag "${project}"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(project), val(meta)

    output:
    tuple val(project), path("*.ped"), emit: ped
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def out   = new File(project + ".ped")
    def samples = (meta.collect().size() > 1) ? meta.sort{ a, b ->
        // First sort on family_id, then on sample id
        a.family_id <=> b.family_id ?: a.id <=> b.id } : meta
    outfile_text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\\t')
    def samples_list = []
    for(int i = 0; i<samples.size(); i++) {
        sample_name = samples[i].id
        if (!samples_list.contains(sample_name)) {
            outfile_text += "\\n" + [samples[i].family_id, sample_name, samples[i].paternal_id, samples[i].maternal_id, samples[i].sex, samples[i].phenotype].join('\\t')
            samples_list.add(sample_name)
        }
    }
    """
    echo -e "$outfile_text" >${project}.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${project}.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_pedigree_file: 1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
