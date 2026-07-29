include { GATK4_CLEANSAM           } from '../../../modules/nf-core/gatk4/cleansam/main'
include { SAMTOOLS_INDEX           } from '../../../modules/nf-core/samtools/index/main'
include { PORTELLO as RUN_PORTELLO } from '../../../modules/nf-core/portello/main'
include { SAMTOOLS_SORT            } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_ADDREPLACERG    } from '../../../modules/nf-core/samtools/addreplacerg/main'
include { SAMTOOLS_CALMD           } from '../../../modules/nf-core/samtools/calmd/main'


workflow PORTELLO {
    take:
    ch_reads_to_assembly_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_assembly_to_ref_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai // channel: [mandatory] [ val(meta), path(fai) ]

    main:

    RUN_PORTELLO(
        ch_assembly_to_ref_bam_bai.join(ch_reads_to_assembly_bam_bai, failOnMismatch: true, failOnDuplicate: true).combine(ch_fasta.map { _meta, fasta -> fasta }).map { meta, asm_to_ref_bam, asm_to_ref_bai, read_to_asm_bam, read_to_asm_bai, ref_fasta ->
            [
                meta,
                asm_to_ref_bam,
                asm_to_ref_bai,
                read_to_asm_bam,
                read_to_asm_bai,
                ref_fasta,
                'partially-phased',
                false,
            ]
        }
    )

    SAMTOOLS_SORT(
        RUN_PORTELLO.out.bam,
        ch_fasta.join(ch_fai).collect(),
        "bai",
    )

    SAMTOOLS_ADDREPLACERG(
        SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.index, failOnMismatch: true, failOnDuplicate: true).map { meta, bam, bai ->
            def read_group = "'@RG\\tID:${meta.id}_hifiasm\\tSM:${meta.id}'"
            [meta, bam, bai, read_group]
        },
        [[], [], [], []],
    )

    // Add MD and NM tags for severus and sniffles
    SAMTOOLS_CALMD(
        SAMTOOLS_ADDREPLACERG.out.bam,
        ch_fasta.join(ch_fai).collect(),
    )

    GATK4_CLEANSAM(
        SAMTOOLS_CALMD.out.bam,
        [[], [], []],
    )

    SAMTOOLS_INDEX(
        GATK4_CLEANSAM.out.bam
    )

    emit:
    bam = GATK4_CLEANSAM.out.bam // channel: [ val(meta), path(bam) ]
    bai = SAMTOOLS_INDEX.out.index // channel: [ val(meta), path(bai) ]
}
