include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'

workflow BAM_TO_FASTQ {

    take:
    ch_sample   // channel: [ val(meta), reads ]

    main:
    ch_versions = Channel.empty()

    // Filter out BAM from fastq
    ch_sample
        .map { meta, fastq -> [ meta + [ 'single_end': true ], fastq ] }
        .branch { meta, reads ->
            fastq: reads.extension == 'gz'
            bam: reads.extension == 'bam'
        }
        .set { ch_filetypes }

    ch_filetypes.fastq.set { ch_sample }

    SAMTOOLS_FASTQ ( ch_filetypes.bam, false )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

    // Mix converted BAM back in
    ch_sample = ch_sample.mix(SAMTOOLS_FASTQ.out.other)

    emit:
    fastq    = ch_sample   // channel: [ val(meta), fastq ]
    versions = ch_versions // channel: [ versions.yml ]
}

