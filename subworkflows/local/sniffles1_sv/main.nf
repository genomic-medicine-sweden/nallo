include { SNIFFLES1 as RUN_SNIFFLES1 } from '../../../modules/local/sniffles1/main'

workflow SNIFFLES1_SV {
    take:
    ch_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]

    main:
    RUN_SNIFFLES1(
        ch_bam_bai
    )

    emit:
    vcf = RUN_SNIFFLES1.out.vcf.map { meta, vcf -> [meta + [sv_caller: 'sniffles1', needs_reheader: true, skip_vep_prep: false], vcf] } // channel: [ val(meta), path(vcf) ]
}
