include { KIVVI as RUN_KIVVI } from '../../../modules/nf-core/kivvi/main'
include { SAMTOOLS_CONVERT   } from '../../../modules/nf-core/samtools/convert/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_MERGE     } from '../../../modules/nf-core/bcftools/merge/main'

workflow KIVVI {
    take:
    ch_bam_bai       // channel: [ val(meta), path(bam), path(bai) ]
    val_kivvi_preset // channel: [ val(meta), path(bed) ]
    val_cram_output  // channel: [ val(meta), path(vcf) ]
    ch_fasta         // channel: [ val(meta), path(fasta) ]
    ch_fai           // channel: [ val(meta), path(fai) ]

    main:
        RUN_KIVVI(
            ch_bam_bai,
            val_kivvi_preset,
        )

        BCFTOOLS_SORT(
            RUN_KIVVI.out.vcf,
        )

        ch_vcfs_per_family = BCFTOOLS_SORT.out.vcf
            .join(BCFTOOLS_SORT.out.tbi, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, vcf, tbi -> [meta.family_id, vcf, tbi] }
            .groupTuple()
            .map { family_id, vcfs, tbis -> [[id: family_id], vcfs, tbis, []] }

        BCFTOOLS_MERGE(
            ch_vcfs_per_family,
            [[], [], []],
        )

        if (val_cram_output) {
            SAMTOOLS_CONVERT(
                RUN_KIVVI.out.bam.join(RUN_KIVVI.out.bai, failOnDuplicate: true, failOnMismatch: true),
                ch_fasta.join(ch_fai, failOnDuplicate: true, failOnMismatch: true).collect(),
            )
        }

    emit:
    bam  = RUN_KIVVI.out.bam                                                                                                 // channel: [ val(meta), path(bam) ]
    bai  = RUN_KIVVI.out.bai                                                                                                 // channel: [ val(meta), path(bai) ]
    crai = val_cram_output ? SAMTOOLS_CONVERT.out.crai : channel.empty()
    cram = val_cram_output ? SAMTOOLS_CONVERT.out.cram : channel.empty()
    svg  = RUN_KIVVI.out.svg                                                                                                 // channel: [ val(meta), path(svg) ]
    json = RUN_KIVVI.out.json
    sample_vcf  = RUN_KIVVI.out.vcf
    sample_tbi  = BCFTOOLS_SORT.out.tbi                                                                                                 // channel: [ val(meta), path(vcf) ]
    family_vcf  = BCFTOOLS_MERGE.out.vcf
    family_tbi  = BCFTOOLS_MERGE.out.index
}
