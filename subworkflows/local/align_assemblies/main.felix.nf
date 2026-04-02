include { MINIMAP2_ALIGN   } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX   } from '../../../modules/nf-core/minimap2/index/main'
include { SAMTOOLS_MERGE   } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_VIEW    } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'
include { TAGBAM           } from '../../../modules/nf-core/tagbam/main'

workflow ALIGN_ASSEMBLIES {

    take:
    ch_assembly // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai      // channel: [mandatory] [ val(meta), path(fai)   ]
    cram_output // bool: Publish alignments as CRAM (true) or BAM (false)

    main:

    MINIMAP2_INDEX (
        ch_fasta
    )

    MINIMAP2_ALIGN (
        ch_assembly.collect().map { meta1, bam1, _meta2, bam2 -> [meta1 - meta1.subMap('haplotype'), [ bam1, bam2] ] },
        MINIMAP2_INDEX.out.index.collect(),
        true,
        'bai',
        false,
        false
    )

    SAMTOOLS_VIEW (
        MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index, failOnMismatch:true, failOnDuplicate:true),
        [[],[]],
        [],
        false
    )

    TAGBAM (
        SAMTOOLS_VIEW.out.bam
    )

    TAGBAM.out.bam
    SAMTOOLS_VIEW.out.bam
        .map { meta, bam -> [ meta - meta.subMap('haplotype'), bam ] }
        .groupTuple(size: 2)
        .set { ch_assemblies_per_sample }

    SAMTOOLS_MERGE (
        ch_assemblies_per_sample,
        [[],[],[],[]],
    )

    // Publish alignment as CRAM if requested
    if (cram_output) {
        SAMTOOLS_CONVERT(
            SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_MERGE.out.index, failOnDuplicate: true, failOnMismatch: true),
            ch_fasta.join(ch_fai, failOnDuplicate: true, failOnMismatch: true).collect(),
        )
    }

    emit:
    bam      = MINIMAP2_ALIGN.out.bam // channel: [ val(meta), path(bam) ]
    bai      = MINIMAP2_ALIGN.out.index // channel: [ val(meta), path(bai) ]
}
