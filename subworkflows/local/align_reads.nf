include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_UNSPLIT        } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SPLIT          } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MINIMAP2_ALIGN } from '../../modules/nf-core/samtools/index/main'
include { FASTP                                           } from '../../modules/nf-core/fastp/main'
include { SAMTOOLS_CAT_SORT_INDEX                         } from '../../modules/local/samtools_cat_sort_index'

workflow ALIGN_READS {

    // Maybe it's possible to do the preprocessing in a separate workflow,
    // then specify in meta if the read that should be aligned is split or not
    // for a cleaner workflow? - branch channel on meta."split"?

    take:
    ch_sample
    ch_mmi   // channel: [ val(meta), mmi ]

    main:
    ch_versions = Channel.empty()
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    ch_bam_bai = Channel.empty()

    // Remap index
    mmi = ch_mmi.map{ it [1] }

    // Split FASTQ
    if (params.split_fastq >= 250) {

        // Add meta info for fastp
        ch_sample
            .map{ meta, fastq -> [ meta + ["single_end":true], fastq]}
            .set{ ch_fastp_in }

        // To run this params.split_fastq should be >= 250
        FASTP( ch_fastp_in, [], [], [] )

        // Transpose and remove single_end from meta - how to just remove one element?
        FASTP.out.reads
            .transpose()
            .map{ meta, split_fastq -> [ [
                'id':meta['id'],
                'family_id':meta['family_id'],
                'paternal_id':meta['paternal_id'],
                'maternal_id':meta['maternal_id'],
                'sex':meta['sex'],
                'phenotype':meta['phenotype'],
                ], split_fastq ]}
            .set { ch_reads }

        MINIMAP2_ALIGN_SPLIT( ch_reads.combine(mmi), true, false, false )

        MINIMAP2_ALIGN_SPLIT.out.bam
            .groupTuple() // Collect aligned files per sample
            .set{ ch_samtools_cat_in }

        // Make one BAM per sample
        SAMTOOLS_CAT_SORT_INDEX(ch_samtools_cat_in)

        // Gather files
        ch_bam = ch_bam.mix(SAMTOOLS_CAT_SORT_INDEX.out.bam)
        ch_bai = ch_bai.mix(SAMTOOLS_CAT_SORT_INDEX.out.bai)
        ch_bam_bai = ch_bam_bai.mix(SAMTOOLS_CAT_SORT_INDEX.out.bam_bai)

        // Gather versions
        ch_versions = ch_versions.mix(FASTP.out.versions)
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_SPLIT.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_CAT_SORT_INDEX.out.versions)

        } else {

        MINIMAP2_ALIGN_UNSPLIT ( ch_sample.combine( mmi ), true, false, false )
        SAMTOOLS_INDEX_MINIMAP2_ALIGN ( MINIMAP2_ALIGN_UNSPLIT.out.bam )

        MINIMAP2_ALIGN_UNSPLIT.out.bam
            .join(SAMTOOLS_INDEX_MINIMAP2_ALIGN.out.bai)
            .set{ ch_bam_bai_single }

        // Gather files
        ch_bam = ch_bam.mix(MINIMAP2_ALIGN_UNSPLIT.out.bam)
        ch_bai = ch_bai.mix(SAMTOOLS_INDEX_MINIMAP2_ALIGN.out.bai)
        ch_bam_bai = ch_bam_bai.mix(ch_bam_bai_single)

        // Gather versions
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_UNSPLIT.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MINIMAP2_ALIGN.out.versions)
    }


    emit:
    bam      = ch_bam       // channel: [ val(meta), bam ]
    bai      = ch_bai       // channel: [ val(meta), bai ]
    bam_bai  = ch_bam_bai   // channel: [ val(meta), bam, bai]
    versions = ch_versions  // channel: [ versions.yml ]
}

