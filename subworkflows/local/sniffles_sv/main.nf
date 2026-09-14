include { SNIFFLES as RUN_SNIFFLES } from '../../../modules/nf-core/sniffles/main'

workflow SNIFFLES_SV {
    take:
    ch_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_tandem_repeats // channel:  [optional] [ val(meta), path(bed) ]

    main:
    RUN_SNIFFLES(
        ch_bam_bai,
        ch_fasta,
        ch_tandem_repeats,
        true,
        false,
    )

    emit:
    vcf = RUN_SNIFFLES.out.vcf.map { meta, vcf -> [meta + [sv_caller: 'sniffles', needs_reheader: true, skip_vep_prep: false], vcf] } // channel: [ val(meta), path(vcf) ]
}
