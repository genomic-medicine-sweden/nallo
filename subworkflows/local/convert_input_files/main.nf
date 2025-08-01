include { SAMTOOLS_IMPORT } from '../../../modules/nf-core/samtools/import/main'
include { SAMTOOLS_FASTQ  } from '../../../modules/nf-core/samtools/fastq/main'

// This subworkflow converts input files between FASTQ and BAM formats,
// while also outputting any input files that were not converted
// into the respective output channels.
workflow CONVERT_INPUT_FILES {
    take:
    ch_input      // channel: [ val(meta), reads ]
    convert_bam   //    bool: Should BAM files be converted to FASTQ
    convert_fastq //    bool: Should FASTQ files be converted to BAM

    main:
    ch_versions = Channel.empty()

    ch_input
        .branch { _meta, reads ->
            fastq: reads.extension == 'gz'
            bam: reads.extension == 'bam'
        }
        .set { reads_to_convert }

    ch_bam = reads_to_convert.bam
    ch_fastq = reads_to_convert.fastq

    if (convert_bam) {
        SAMTOOLS_FASTQ(
            reads_to_convert.bam,
            false,
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

        // Mix converted files back in
        ch_fastq = ch_fastq.mix(SAMTOOLS_FASTQ.out.other)
    }
    if (convert_fastq) {
        SAMTOOLS_IMPORT(
            reads_to_convert.fastq
        )
        ch_versions = ch_versions.mix(SAMTOOLS_IMPORT.out.versions)

        // Mix converted files back in
        ch_bam = ch_bam.mix(SAMTOOLS_IMPORT.out.bam)
    }

    emit:
    bam      = ch_bam // channel: [ val(meta), bam ]
    fastq    = ch_fastq // channel: [ val(meta), fastq ]
    versions = ch_versions // channel: [ versions.yml ]
}
