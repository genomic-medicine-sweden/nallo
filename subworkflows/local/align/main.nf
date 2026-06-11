include { MINIMAP2_ALIGN   } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE   } from '../../../modules/nf-core/samtools/merge/main'
include { SPLITUBAM        } from '../../../modules/nf-core/splitubam/main'

workflow ALIGN {

    take:
    ch_ubam             // channel: [mandatory] [ val(meta), path(reads)]
    ch_mmi              // channel: [mandatory] [ val(meta), path(index) ]
    val_split_alignment // bool

    main:

    if (val_split_alignment) {
        SPLITUBAM(ch_ubam)
        ch_unmapped = SPLITUBAM.out.bam
            .transpose()
    } else {
        ch_unmapped = ch_ubam
    }

    // Adding file key to meta for proper joining of minimap output
    ch_minimap_in = ch_unmapped
        .map { meta, bam -> tuple(meta + [file: bam.name], bam) }

    /*
     * Create a grouping key per sample that records the number of split files,
     * allowing downstream merging to trigger as soon as all alignments of a sample are ready.
     */
    ch_reads_grouping_key = ch_unmapped
        .groupTuple()
        .map { meta, files -> tuple(meta.id, files.size()) }

    MINIMAP2_ALIGN(
        ch_minimap_in,
        ch_mmi,
        true,
        'bai',
        false,
        false,
    )

    ch_mapped = MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index, failOnMismatch: true, failOnDuplicate:true)
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

    SAMTOOLS_MERGE(
        ch_mapped,
        [[], [], [], []]
    )

    emit:
    bam = SAMTOOLS_MERGE.out.bam.map { meta, bam -> tuple(meta - meta.subMap('n_files'), bam) }
    bai = SAMTOOLS_MERGE.out.index.map { meta, bai -> tuple(meta - meta.subMap('n_files'), bai) }
}
