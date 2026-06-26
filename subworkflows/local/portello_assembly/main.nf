include { GATK4_CLEANSAM                            } from '../../../modules/nf-core/gatk4/cleansam/main'
include { PBMM2_ALIGN                               } from '../../../modules/nf-core/pbmm2/align/main'
include { FIND_CONCATENATE                          } from '../../../modules/nf-core/find/concatenate/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PBMM2    } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CLEANSAM } from '../../../modules/nf-core/samtools/index/main'
include { PORTELLO                                  } from '../../../modules/nf-core/portello/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_PORTELLO   } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_ADDREPLACERG                     } from '../../../modules/nf-core/samtools/addreplacerg/main'
include { SAMTOOLS_CALMD                            } from '../../../modules/nf-core/samtools/calmd/main'


workflow PORTELLO_ASSEMBLY {
    take:
    ch_assembled_haplotypes     // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_bam                      // channel: [mandatory] [ val(meta), path(bam) ]
    ch_assembly_bam_bai         // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta                    // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai                      // channel: [mandatory] [ val(meta), path(fai) ]

    main:
    ch_assembled_haplotypes
        .map { meta, bam -> [ meta - meta.subMap('haplotype'), bam ] }
        .groupTuple(size: 2)
        .set { ch_assemblies_to_concatenate }

    FIND_CONCATENATE (
        ch_assemblies_to_concatenate
    )

    // Match BAM files and reference by sample ID for alignment with pbmm2
    ch_bam
        .map { meta, bam -> [[id: meta.id], meta, bam] }
        .join(
            FIND_CONCATENATE.out.file_out.map { meta, ref -> [[id: meta.id], ref] },
            failOnMismatch: true,
            failOnDuplicate: true
        )
        .multiMap { _key, meta, bam, ref ->
            reads:     [meta, bam]
            reference: [meta, ref]
        }
        .set { ch_pbmm2_input }

    PBMM2_ALIGN(
        ch_pbmm2_input.reads,
        ch_pbmm2_input.reference,
    )

    SAMTOOLS_INDEX_PBMM2(PBMM2_ALIGN.out.bam)

    PORTELLO (
        ch_assembly_bam_bai
            .join(PBMM2_ALIGN.out.bam)
            .join(SAMTOOLS_INDEX_PBMM2.out.index)
            .combine(ch_fasta.map { _meta, fasta -> fasta })
            .map { meta, asm_to_ref_bam, asm_to_ref_bai, read_to_asm_bam, read_to_asm_bai, ref_fasta ->
                [
                    meta,
                    asm_to_ref_bam,
                    asm_to_ref_bai,
                    read_to_asm_bam,
                    read_to_asm_bai,
                    ref_fasta,
                    'partially-phased',
                    false
                ]
            }
    )

    SAMTOOLS_SORT_PORTELLO (
        PORTELLO.out.bam,
        ch_fasta.join(ch_fai).collect(),
        "bai"
    )

    SAMTOOLS_ADDREPLACERG (
        SAMTOOLS_SORT_PORTELLO.out.bam
            .join(SAMTOOLS_SORT_PORTELLO.out.index, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, bam, bai ->
                def read_group = "'@RG\\tID:${meta.id}_hifiasm\\tSM:${meta.id}'"
                [meta, bam, bai, read_group]
            },
        [[],[],[],[]]
    )

    // Add MD and NM tags for severus
    SAMTOOLS_CALMD (
        SAMTOOLS_ADDREPLACERG.out.bam,
        ch_fasta.join(ch_fai).collect(),
    )

    GATK4_CLEANSAM(
        SAMTOOLS_CALMD.out.bam,
        [[],[],[]],
    )

    SAMTOOLS_INDEX_CLEANSAM(
        GATK4_CLEANSAM.out.bam
    )

    emit:
    bam = GATK4_CLEANSAM.out.bam            // channel: [ val(meta), path(bam) ]
    bai = SAMTOOLS_INDEX_CLEANSAM.out.index // channel: [ val(meta), path(bai) ]
}
