include { MINIMAP2_ALIGN   } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE   } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_VIEW    } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'
include { TAGBAM           } from '../../../modules/nf-core/tagbam/main'

workflow ALIGN_ASSEMBLIES {
    take:
    ch_assembly // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_mmi      // channel: [mandatory] [ val(meta), path(mmi) ]

    main:

    MINIMAP2_ALIGN(
        ch_assembly,
        ch_mmi,
        true,
        'bai',
        false,
        false,
    )

    SAMTOOLS_VIEW(
        MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index, failOnMismatch: true, failOnDuplicate: true),
        [[], [], []],
        [],
        false,
    )

    TAGBAM(
        SAMTOOLS_VIEW.out.bam
    )

    ch_assemblies_per_sample = TAGBAM.out.bam
        .map { meta, bam -> [meta - meta.subMap('haplotype'), bam] }
        .groupTuple(size: 2)
        .map { meta, bams -> [meta, bams, []] }

    SAMTOOLS_MERGE(
        ch_assemblies_per_sample,
        [[], [], [], []],
    )

    emit:
    bam  = SAMTOOLS_MERGE.out.bam                                    // channel: [ val(meta), path(bam) ]
    bai  = SAMTOOLS_MERGE.out.index                                  // channel: [ val(meta), path(bai) ]
    cram = cram_output ? SAMTOOLS_CONVERT.out.cram : channel.empty() // channel: [ val(meta), path(cram) ]
    crai = cram_output ? SAMTOOLS_CONVERT.out.crai : channel.empty() // channel: [ val(meta), path(crai) ]
}
