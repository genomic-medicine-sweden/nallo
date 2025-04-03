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
    ch_versions = Channel.empty()

    MINIMAP2_INDEX (
        ch_fasta
    )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    MINIMAP2_ALIGN (
        ch_assembly,
        MINIMAP2_INDEX.out.index.collect(),
        true,
        'bai',
        false,
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    SAMTOOLS_VIEW (
        MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index, failOnMismatch:true, failOnDuplicate:true),
        [[],[]],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    TAGBAM (
        SAMTOOLS_VIEW.out.bam
    )
    ch_versions = ch_versions.mix(TAGBAM.out.versions)

    TAGBAM.out.bam
        .map { meta, bam -> [ meta - meta.subMap('haplotype'), bam ] }
        .groupTuple(size: 2)
        .set { ch_assemblies_per_sample }

    SAMTOOLS_MERGE (
        ch_assemblies_per_sample,
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // Publish alignment as CRAM if requested
    if (cram_output) {
        SAMTOOLS_CONVERT(
            SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_MERGE.out.bai, failOnDuplicate: true, failOnMismatch: true),
            ch_fasta,
            ch_fai
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)
    }

    emit:
    bam      = SAMTOOLS_MERGE.out.bam // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_MERGE.out.bai // channel: [ val(meta), path(bai) ]
    versions = ch_versions            // channel: [ versions.yml ]
}

