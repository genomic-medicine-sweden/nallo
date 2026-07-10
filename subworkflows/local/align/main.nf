include { MINIMAP2_ALIGN   } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX   } from '../../../modules/nf-core/minimap2/index/main'

workflow ALIGN {

    take:
    ch_ubam             // channel: [mandatory] [ val(meta), path(reads)]
    ch_fasta            // channel: [mandatory] [ val(meta), path(fasta) ]

    main:

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

    emit:
    bam   = MINIMAP2_ALIGN.out.bam
    index = MINIMAP2_ALIGN.out.index
}
