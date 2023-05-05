//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow PED_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:['family_id', 'id', 'paternal_id', 'maternal_id', 'sex', 'phenotype'], sep:'\t' )
        .map { create_ped_channel(it) }
        .set { ch_ped_processed }

    emit:
    ch_ped_processed                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_ped_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id                = row.id
    meta.family_id                = row.family_id
    meta.paternal_id                = row.paternal_id
    meta.maternal_id                = row.maternal_id
    meta.sex                = row.sex
    meta.phenotype                = row.phenotype
    // add path(s) of the fastq file(s) to the meta map
    def ped_meta = []
    
    ped_meta = meta
    
    return ped_meta
}
