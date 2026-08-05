include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../../modules/nf-core/minimap2/index/main'
include { MM2PLUS_ALIGN  } from '../../../modules/nf-core/mm2plus/align/main'
include { MM2PLUS_INDEX  } from '../../../modules/nf-core/mm2plus/index/main'
include { PBMM2_ALIGN    } from '../../../modules/nf-core/pbmm2/align/main'

workflow ALIGN {
    take:
    ch_ubam // channel: [mandatory] [ val(meta), path(reads)]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    val_aligner // channel: [mandatory] [ val(aligner) ]
    use_genome_reference // boolean: [mandatory] [ true|false ]

    main:

    if (val_aligner == 'pbmm2') {
        // Match BAM files and reference by sample ID if needed (for portello)
        ch_pbmm2_input = ch_ubam
            .combine(ch_fasta)
            .filter { bam_meta, _bam, ref_meta, _ref -> use_genome_reference || bam_meta.id == ref_meta.id }
            .multiMap { bam_meta, bam, _ref_meta, ref ->
                reads: [bam_meta, bam]
                reference: [bam_meta, ref]
            }

        PBMM2_ALIGN(
            ch_pbmm2_input.reads,
            ch_pbmm2_input.reference,
        )

        ch_aligned_reads_bam = PBMM2_ALIGN.out.bam
        ch_aligned_reads_bai = PBMM2_ALIGN.out.index
    }
    if (val_aligner == 'mm2plus') {

        MM2PLUS_INDEX(
            ch_fasta
        )

        MM2PLUS_ALIGN(
            ch_ubam,
            MM2PLUS_INDEX.out.index.collect(),
            true,
            'bai',
            false,
            false,
        )

        ch_aligned_reads_bam = MM2PLUS_ALIGN.out.bam
        ch_aligned_reads_bai = MM2PLUS_ALIGN.out.index
    }
    if (val_aligner == 'minimap2') {

        MINIMAP2_INDEX(
            ch_fasta
        )

        MINIMAP2_ALIGN(
            ch_ubam,
            MINIMAP2_INDEX.out.index.collect(),
            true,
            'bai',
            false,
            false,
        )

        ch_aligned_reads_bam = MINIMAP2_ALIGN.out.bam
        ch_aligned_reads_bai = MINIMAP2_ALIGN.out.index
    }

    emit:
    bam   = ch_aligned_reads_bam
    index = ch_aligned_reads_bai
}
