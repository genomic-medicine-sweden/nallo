include { TRGT                                  } from '../../modules/local/trgt'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRGT } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRGT   } from '../../modules/nf-core/samtools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_TRGT   } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_TRGT } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_MERGE                        } from '../../modules/nf-core/bcftools/merge/main'

workflow REPEAT_ANALYSIS {

    take:
    ch_bam_bai
    ch_fasta
    ch_fai
    ch_trgt_bed

    main:
    ch_repeat_calls_vcf = Channel.empty()
    ch_versions         = Channel.empty()

    ch_bam_bai
        .map{ meta, bam, bai -> [meta, bam, bai, meta.sex] }
        .set{ ch_trgt_input }

    // Run TGRT
    TRGT ( ch_trgt_input, ch_fasta, ch_trgt_bed )

    // Sort and index bam
    SAMTOOLS_SORT_TRGT ( TRGT.out.bam, [[],[]] )
    SAMTOOLS_INDEX_TRGT(SAMTOOLS_SORT_TRGT.out.bam)

    // Sort and index bcf
    BCFTOOLS_SORT_TRGT(TRGT.out.vcf)
    BCFTOOLS_INDEX_TRGT(BCFTOOLS_SORT_TRGT.out.vcf)

    BCFTOOLS_SORT_TRGT.out.vcf
        .join(BCFTOOLS_INDEX_TRGT.out.csi)
        .toList()
        .filter { it.size() > 1 }
        .flatMap()
        .map { meta, bcf, csi -> [ [ id : 'multisample' ], bcf, csi ] }
        .groupTuple()
        .set{ ch_bcftools_merge_in }

    BCFTOOLS_MERGE ( ch_bcftools_merge_in, ch_fasta, ch_fai, [] )

    ch_versions = ch_versions.mix(TRGT.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_TRGT.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TRGT.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_TRGT.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_TRGT.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

