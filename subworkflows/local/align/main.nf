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
        ch_ubam
             .map { meta, bam, _bai -> [meta, bam] }
             .set { ch_splitubam_in }
        SPLITUBAM(ch_splitubam_in)
        ch_unmapped = SPLITUBAM.out.bam
    } else {
        ch_unmapped = ch_ubam
    }

    /*
     * Create a grouping key per sample that records the number of split files,
     * allowing downstream merging to trigger as soon as all alignments of a sample are ready.
     */
    ch_unmapped
        .mix(ch_mapped_bam)
        // Unmapped BAMs have no index but mapped ones do. To avoid an extra processing step, we can just assign a default value to the index.
        .map { meta, bam, _index=null -> [meta - meta.subMap('file'), bam] }
        .groupTuple()
        .map { meta, files -> [meta.id, files.size() ] }
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
            [bam_meta - bam_meta.subMap('file') + [n_files: group_size], bam, bai]
        }
        .map { meta, bam, bai ->
            [groupKey(meta, meta.n_files), bam, bai]
        }
        .groupTuple()
        .map { key, bams, bais -> [key.getGroupTarget(), bams, bais ]}
        .set { ch_mapped }

    SAMTOOLS_MERGE(
        ch_mapped,
        [[], [], [], []]
    )

    SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_MERGE.out.index, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, bam, bai -> [meta - meta.subMap('n_files'), bam, bai] }
        .set { ch_aligned_bam }


    emit:
    bam_bai = ch_aligned_bam
}
