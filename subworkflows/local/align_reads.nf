// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main.nf'

workflow ALIGN_READS {

    take:
    // TODO nf-core: edit input (take) channels
    ch_fastq // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()
    
    MINIMAP2_INDEX( ch_fasta )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())
    
    MINIMAP2_ALIGN( ch_fastq.combine(MINIMAP2_INDEX.out.index.map{ it[1] }), true, false, false)

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    SAMTOOLS_INDEX ( MINIMAP2_ALIGN.out.bam )

    MINIMAP2_ALIGN.out.bam
    .concat(SAMTOOLS_INDEX.out.bai)
    .groupTuple().flatten().collate(3)
    .set{ch_bam_bai}

    emit:
    bam_bai = ch_bam_bai // channel: [ [meta], bam, bai]
    versions = ch_versions                     // channel: [ versions.yml ]
}

