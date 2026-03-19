process CREATE_PEDIGREE_FILE {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), val(sample_metas)

    output:
    tuple val(meta), path("*.ped"), emit: ped
    tuple val("${task.process}"), val('create_pedigree_file'), val("1.0"), emit: versions_create_pedigree_file, topic: versions
    tuple val("${task.process}"), val('create_pedigree_file'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def samples = (sample_metas.collect().size() > 1) ? sample_metas.sort{ a, b ->
        // First sort on family_id, then on sample id
        a.family_id <=> b.family_id ?: a.id <=> b.id } : sample_metas
    outfile_text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\\t')
    def samples_list = []
    samples.each { sample ->
        if (!samples_list.contains(sample.id)) {
            outfile_text += "\\n" + [sample.family_id, sample.id, sample.paternal_id, sample.maternal_id, sample.sex, sample.phenotype].join('\\t')
            samples_list.add(sample.id)
        }
    }
    """
    echo -e "$outfile_text" > ${prefix}.ped
    """

    stub:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ped
    """
}
