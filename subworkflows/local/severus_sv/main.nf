include { SEVERUS as RUN_SEVERUS } from '../../../modules/nf-core/severus/main'

workflow SEVERUS_SV {
    take:
    ch_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_tandem_repeats // channel:  [optional] [ val(meta), path(bed) ]

    main:
    RUN_SEVERUS(
        ch_bam_bai.map { meta, bam, bai -> [meta, bam, bai, [], [], []] },
        ch_tandem_repeats,
    )

    emit:
    vcf = RUN_SEVERUS.out.all_vcf.map { meta, vcf -> [meta + [sv_caller: 'severus', needs_reheader: true, skip_vep_prep: false], vcf] } // channel: [ val(meta), path(vcf) ]
}
