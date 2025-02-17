include { DIPCALL                                  } from '../../../modules/local/dipcall'
include { MINIMAP2_ALIGN                           } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX                           } from '../../../modules/nf-core/minimap2/index/main'
include { SAMTOOLS_INDEX                           } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                           } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_VIEW                            } from '../../../modules/nf-core/samtools/view/main'
include { TAGBAM                                   } from '../../../modules/nf-core/tagbam/main'

workflow ASSEMBLY_VARIANT_CALLING {

    take:
    ch_haplotypes     // channel: [mandatory] [ val(meta), path(paternal_haplotype), path(maternal_haplotype) ]
    ch_paternal_fasta // channel: [mandatory] [ val(meta), path(paternal_fasta) ]
    ch_maternal_fasta // channel: [mandatory] [ val(meta), path(maternal_fasta) ]
    ch_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai            // channel: [mandatory] [ val(meta), path(fai) ]
    ch_par            // channel: |mandatory] [ val(meta), path(par) ]

    main:
    ch_sv_calls_vcf = Channel.empty()
    ch_versions     = Channel.empty()

    ch_haplotypes
        .map{ meta, hap1, hap2 -> [meta, hap1, hap2, meta.sex] }
        .set{ ch_dipcall_input }

    MINIMAP2_INDEX (
        ch_fasta
    )

    MINIMAP2_ALIGN (
        ch_paternal_fasta.mix(ch_maternal_fasta),
        MINIMAP2_INDEX.out.index,
        true,
        'bai',
        false,
        false
    )

    SAMTOOLS_VIEW (
        MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index),
        [[],[]],
        []
    )

    TAGBAM (
        SAMTOOLS_VIEW.out.bam
    )

    TAGBAM.out.bam
        .map { meta, bam -> [ [ meta - meta.subMap('haplotype') ], bam ] }
        .groupTuple()
        .view()

    SAMTOOLS_MERGE (
        TAGBAM.out.bam,
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    ch_sv_calls_vcf        // channel: ?
    versions = ch_versions // channel: [ versions.yml ]
}

