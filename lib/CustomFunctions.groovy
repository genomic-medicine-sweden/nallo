import nextflow.Nextflow

class CustomFunctions {

    // Function to generate a pedigree file
    public static File makePed(samples, outdir) {
        def case_name  = "multisample"
        def outfile  = new File(outdir +"/pipeline_info/${case_name}" + '.ped')
        outfile.text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\t')
        def samples_list = []
        for(int i = 0; i<samples.size(); i++) {
            samples[i] = samples[i][0]
            def sample_name =  samples[i].id
            if (!samples_list.contains(sample_name)) {
                outfile.append('\n' + [samples[i].family_id, sample_name, samples[i].paternal_id, samples[i].maternal_id, samples[i].sex, samples[i].phenotype].join('\t'));
                samples_list.add(sample_name)
            }
        }
        return outfile
    }
}
