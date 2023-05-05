//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { ch_sample }

    emit:
    ch_sample                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.file).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> fastq.gz file does not exist!\n${row.file}"
    }
    
    fastq_meta = [ meta, file(row.file) ]
    
    return fastq_meta
}
