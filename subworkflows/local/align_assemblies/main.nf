include { MINIMAP2_ALIGN   } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX   } from '../../../modules/nf-core/minimap2/index/main'
include { MM2PLUS_ALIGN    } from '../../../modules/nf-core/mm2plus/align/main'
include { SAMTOOLS_MERGE   } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_VIEW    } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'
include { TAGBAM           } from '../../../modules/nf-core/tagbam/main'

workflow ALIGN_ASSEMBLIES {
    take:
    ch_assembly // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai // channel: [mandatory] [ val(meta), path(fai)   ]
    val_cram_output // bool: Publish alignments as CRAM (true) or BAM (false)
    val_assembly_aligner // string: Which aligner to use for assembly alignment

    main:

    MINIMAP2_INDEX(
        ch_fasta
    )

    if (val_assembly_aligner == 'mm2plus') {

        MM2PLUS_ALIGN(
            ch_assembly,
            MM2PLUS_INDEX.out.index.collect(),
            true,
            'bai',
            false,
            false,
        )
        ch_aligned_assemblies_bam = MM2PLUS_ALIGN.out.bam
    }
    else {

        MINIMAP2_ALIGN(
            ch_assembly,
            MINIMAP2_INDEX.out.index.collect(),
            true,
            'bai',
            false,
            false,
        )

        ch_aligned_assemblies_bam = MINIMAP2_ALIGN.out.bam
    }

    TAGBAM(
        ch_aligned_assemblies_bam
    )

    ch_assemblies_per_sample = TAGBAM.out.bam
        .map { meta, bam -> [meta - meta.subMap('haplotype'), bam] }
        .groupTuple(size: 2)
        .map { meta, bams -> [meta, bams, []] }

    SAMTOOLS_MERGE(
        ch_assemblies_per_sample,
        [[], [], [], []],
    )

    SAMTOOLS_VIEW(
        SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_MERGE.out.index, failOnMismatch: true, failOnDuplicate: true),
        [[], [], []],
        [],
        'bai',
    )

    // Publish alignment as CRAM if requested
    if (val_cram_output) {
        SAMTOOLS_CONVERT(
            SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_VIEW.out.bai, failOnDuplicate: true, failOnMismatch: true),
            ch_fasta.join(ch_fai, failOnDuplicate: true, failOnMismatch: true).collect(),
        )
    }

    emit:
    unfiltered_bam = SAMTOOLS_MERGE.out.bam // channel: [ val(meta), path(bam) ]
    unfiltered_bai = SAMTOOLS_MERGE.out.index // channel: [ val(meta), path(bai) ]
    bam            = SAMTOOLS_VIEW.out.bam // channel: [ val(meta), path(bam) ]
    bai            = SAMTOOLS_VIEW.out.bai // channel: [ val(meta), path(bai) ]
    cram           = val_cram_output ? SAMTOOLS_CONVERT.out.cram : channel.empty() // channel: [ val(meta), path(cram) ]
    crai           = val_cram_output ? SAMTOOLS_CONVERT.out.crai : channel.empty() // channel: [ val(meta), path(crai) ]
}
