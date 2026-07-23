include { MINIMAP2_ALIGN } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../../modules/nf-core/minimap2/index/main'
include { PBMM2_ALIGN    } from '../../../modules/nf-core/pbmm2/align/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_CALMD } from '../../../modules/nf-core/samtools/calmd/main'

workflow ALIGN {
    take:
    ch_ubam // channel: [mandatory] [ val(meta), path(reads)]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai // channel: [mandatory] [ val(meta), path(fai) ]
    val_aligner // channel: [mandatory] [ val(aligner) ]
    val_skip_portello // boolean: [mandatory] [ true|false ]

    main:

    if (val_aligner == 'pbmm2') {

        if (!val_skip_portello) {
            // Match BAM files and reference by sample ID for alignment with pbmm2
            ch_ubam
                .map { meta, bam -> [[id: meta.id], meta, bam] }
                .combine(
                    ch_fasta.map { meta, ref -> [[id: meta.id], ref] },
                    by: 0
                )
                .multiMap { _id, meta, bam, ref ->
                    reads: [meta, bam]
                    reference: [meta, ref]
                }
                .set { ch_pbmm2_input }
        }
        else {
            ch_pbmm2_input = ch_ubam
                .combine(ch_fasta)
                .multiMap { bam_meta, bam, _ref_meta, ref ->
                    reads: [bam_meta, bam]
                    reference: [bam_meta, ref]
                }
        }

        PBMM2_ALIGN(
            ch_pbmm2_input.reads,
            ch_pbmm2_input.reference,
        )

        if (val_skip_portello) {
            // Add MD and NM tags for sniffles and severus
            SAMTOOLS_CALMD(
                PBMM2_ALIGN.out.bam,
                ch_fasta.join(ch_fai).collect(),
            )

            ch_aligned_reads_bam = SAMTOOLS_CALMD.out.bam
        }
        else {
            ch_aligned_reads_bam = PBMM2_ALIGN.out.bam
        }

        SAMTOOLS_INDEX(ch_aligned_reads_bam)

        ch_aligned_reads_bai = SAMTOOLS_INDEX.out.index
    }
    else {

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
