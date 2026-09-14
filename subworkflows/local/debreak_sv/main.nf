include { DEBREAK as RUN_DEBREAK } from '../../../modules/nf-core/debreak/main'

workflow DEBREAK_SV {
    take:
    ch_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]

    main:
    RUN_DEBREAK(
        ch_bam_bai,
        ch_fasta,
    )

    emit:
    vcf = RUN_DEBREAK.out.vcf.map { meta, vcf -> [meta + [sv_caller: 'debreak', needs_reheader: true, skip_vep_prep: false], vcf] } // channel: [ val(meta), path(vcf) ]
}
