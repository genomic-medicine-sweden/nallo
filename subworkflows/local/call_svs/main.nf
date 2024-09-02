include { SNIFFLES as SNIFFLES_MULTISAMPLE } from '../../../modules/nf-core/sniffles/main'
include { SNIFFLES                         } from '../../../modules/nf-core/sniffles/main'

workflow CALL_SVS {

    take:
    ch_bam_bai        // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai            // channel: [mandatory] [ val(meta), path(fai) ]
    ch_tandem_repeats // channel: [optional]  [ val(meta), path(bed) ]

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
    ch_sv_calls_vcf = SNIFFLES.out.vcf             // channel: [ val(meta), path(vcf) ]
    ch_multisample  = SNIFFLES_MULTISAMPLE.out.vcf // channel: [ val(meta), path(vcf) ]
    versions        = ch_versions                  // channel: [ path(versions.yml) ]
}

