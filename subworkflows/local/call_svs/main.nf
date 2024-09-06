include { BCFTOOLS_MERGE                         } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_REHEADER as REHEADER_SNIFFLES } from '../../../modules/nf-core/bcftools/reheader/main'
include { CREATE_SAMPLES_FILE                    } from '../../../modules/local/create_samples_file.nf'
include { SNIFFLES                               } from '../../../modules/nf-core/sniffles/main'

workflow CALL_SVS {

    take:
    ch_bam_bai        // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai            // channel: [mandatory] [ val(meta), path(fai) ]
    ch_tandem_repeats // channel: [optional]  [ val(meta), path(bed) ]

    main:
    ch_versions     = Channel.empty()

    SNIFFLES (ch_bam_bai, ch_fasta, ch_tandem_repeats, true, false)
    ch_versions = ch_versions.mix(SNIFFLES.out.versions)

    CREATE_SAMPLES_FILE ( SNIFFLES.out.vcf.map { meta, vcf -> meta } )
    ch_versions = ch_versions.mix(CREATE_SAMPLES_FILE.out.versions)

    SNIFFLES.out.vcf
        .join( CREATE_SAMPLES_FILE.out.samples )
        .map { meta, vcf, samples -> [ meta, vcf, [], samples] }
        .set { ch_reheader_sniffles_in }

    REHEADER_SNIFFLES ( ch_reheader_sniffles_in, [[],[]] )
    ch_versions = ch_versions.mix(REHEADER_SNIFFLES.out.versions)

    REHEADER_SNIFFLES.out.vcf
        .join(REHEADER_SNIFFLES.out.index)
        .map { meta, vcf, tbi -> [ [ 'id': meta.project ], vcf, tbi ] }
        .groupTuple()
        .set{ ch_bcftools_merge_in }

    BCFTOOLS_MERGE ( ch_bcftools_merge_in, ch_fasta, ch_fai, [[],[]] )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    emit:
    ch_sv_calls_vcf     = REHEADER_SNIFFLES.out.vcf   // channel: [ val(meta), path(vcf) ]
    ch_sv_calls_tbi     = REHEADER_SNIFFLES.out.index // channel: [ val(meta), path(tbi) ]
    ch_multisample_vcf  = BCFTOOLS_MERGE.out.vcf      // channel: [ val(meta), path(vcf) ]
    ch_multisample_tbi  = BCFTOOLS_MERGE.out.index    // channel: [ val(meta), path(tbi) ]
    versions            = ch_versions                 // channel: [ path(versions.yml) ]
}

