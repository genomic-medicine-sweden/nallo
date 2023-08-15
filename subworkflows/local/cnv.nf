include { HIFICNV } from '../../modules/local/hificnv'

workflow CNV {

    take:
    ch_bam_bai // channel: [ val(meta), [[ bam ], [bai]] ]
    ch_fasta

    main:
    ch_versions     = Channel.empty()

    HIFICNV(ch_bam_bai, ch_fasta, [[],[]], [], [])

    ch_versions = ch_versions.mix(HIFICNV.out.versions)

    emit:
    versions        = ch_versions                  // channel: [ versions.yml ]
}

