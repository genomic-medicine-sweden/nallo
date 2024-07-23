include { SNIFFLES as SNIFFLES_MULTISAMPLE } from '../../modules/nf-core/sniffles/main'
include { SNIFFLES                         } from '../../modules/nf-core/sniffles/main'

workflow STRUCTURAL_VARIANT_CALLING {

    take:
    ch_bam_bai // channel: [ val(meta), [[ bam ], [bai]] ]
    ch_fasta
    ch_fai
    ch_tandem_repeats

    main:
    ch_versions     = Channel.empty()

    SNIFFLES (ch_bam_bai, ch_fasta, ch_tandem_repeats, true, true)

    SNIFFLES.out.snf
        .map { meta, snf -> [ [ 'id': meta.project ], snf ] }
        .groupTuple()
        .map { meta, snfs -> [ meta, snfs, [] ] }
        .set{ ch_multisample_input }

    SNIFFLES_MULTISAMPLE( ch_multisample_input, ch_fasta, ch_tandem_repeats, true, false )

    ch_versions = ch_versions.mix(SNIFFLES.out.versions)
    ch_versions = ch_versions.mix(SNIFFLES_MULTISAMPLE.out.versions)

    emit:
    ch_sv_calls_vcf = SNIFFLES.out.vcf             // channel: [ val(meta), vcf]
    ch_multisample  = SNIFFLES_MULTISAMPLE.out.vcf // channel: [ val(meta), multisample.sniffles.vcf ]
    versions        = ch_versions                  // channel: [ versions.yml ]
}

