include { MINIMAP2_ALIGN   } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE   } from '../../../modules/nf-core/samtools/merge/main'
include { SPLITUBAM        } from '../../../modules/nf-core/splitubam/main'

workflow ALIGN {

    take:
    ch_ubam             // channel: [mandatory] [ val(meta), path(reads)]
    ch_mapped_bam       // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_mmi              // channel: [mandatory] [ val(meta), path(fasta) ]
    val_split_alignment // bool

    main:

    if (val_split_alignment) {
        SPLITUBAM(ch_ubam)
        SPLITUBAM.out.bam
            .transpose()
            // Adding file key to meta for proper joining of minimap output
            .map { meta, bam -> tuple(meta + [file: bam.name], bam) }
            .set { ch_unmapped }
    } else {
        ch_unmapped = ch_ubam
    }

    /*
     * Create a grouping key per sample that records the number of split files,
     * allowing downstream merging to trigger as soon as all alignments of a sample are ready.
     */
    ch_unmapped
        // add dummy index for unmapped BAMs to allow mixing with mapped BAMs
        .map { meta, bam -> tuple(meta - meta.subMap('file'), bam, [])}
        .mix(ch_mapped_bam)
        .groupTuple()
        .map { meta, files, _indexes -> tuple(meta.id, files.size()) }
        .set { ch_reads_grouping_key }

    MINIMAP2_ALIGN(
        ch_unmapped,
        ch_mmi,
        true,
        'bai',
        false,
        false,
    )

    MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index, failOnMismatch: true, failOnDuplicate:true)
        .mix(ch_mapped_bam)
        .combine(ch_reads_grouping_key)
        .filter { bam_meta, _bam, _bai, group_id, _group_size ->
            bam_meta.id == group_id
        }
        .map { bam_meta, bam, bai, _group_id, group_size ->
            tuple(bam_meta - bam_meta.subMap('file') + [n_files: group_size], bam, bai)
        }
        .map { meta, bam, bai ->
            tuple(groupKey(meta, meta.n_files), bam, bai)
        }
        .groupTuple()
        .map { key, bams, bais -> tuple(key.getGroupTarget(), bams, bais) }
        .set { ch_mapped }

    SAMTOOLS_MERGE(
        ch_mapped,
        [[], [], [], []]
    )

    SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_MERGE.out.index, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, bam, bai -> tuple(meta - meta.subMap('n_files'), bam, bai) }
        .set { ch_aligned_bam }


    emit:
    bam_bai = ch_aligned_bam
}
