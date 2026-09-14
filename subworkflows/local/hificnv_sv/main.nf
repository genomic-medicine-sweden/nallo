include { HIFICNV as RUN_HIFICNV       } from '../../../modules/nf-core/hificnv/main'
include { TABIX_TABIX as TABIX_HIFICNV } from '../../../modules/nf-core/tabix/tabix/main'

workflow HIFICNV_SV {
    take:
    ch_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_snvs // channel:  [optional] [ val(meta), path(vcf) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_expected_xy_bed // channel: [mandatory] [ val(meta), path(bed) ]
    ch_expected_xx_bed // channel: [mandatory] [ val(meta), path(bed) ]
    ch_exclude_bed // channel:  [optional] [ val(meta), path(bed) ]
    create_maf_track //     bool: [mandatory] should HiFiCNV produce a MAF track

    main:
    ch_for_hificnv = channel.empty()
    if (create_maf_track) {
        ch_for_hificnv = ch_bam_bai.join(ch_snvs, failOnMismatch: true, failOnDuplicate: true)
    }
    else {
        ch_for_hificnv = ch_bam_bai.map { meta, bam, bai -> [meta, bam, bai, []] }
    }

    ch_hificnv_input = ch_for_hificnv
        .combine(ch_expected_xy_bed)
        .combine(ch_expected_xx_bed)
        .multiMap { meta, bam, bai, maf, xy_meta, xy_bed, xx_meta, xx_bed ->
            def expected_cn_meta = meta.sex == 1 ? xy_meta : xx_meta
            def expected_cn_bed = meta.sex == 1 ? xy_bed : xx_bed
            bam_bai_maf: [meta, bam, bai, maf]
            expected_cn: [expected_cn_meta, expected_cn_bed]
        }

    RUN_HIFICNV(
        ch_hificnv_input.bam_bai_maf,
        ch_fasta,
        ch_exclude_bed,
        ch_hificnv_input.expected_cn,
    )

    TABIX_HIFICNV(
        RUN_HIFICNV.out.vcf
    )

    ch_vcf_tbi = RUN_HIFICNV.out.vcf
        .join(TABIX_HIFICNV.out.index, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi -> [meta + [sv_caller: 'hificnv', needs_reheader: false, skip_vep_prep: true], vcf, tbi] }

    emit:
    vcf     = ch_vcf_tbi.map { meta, vcf, _tbi -> [meta, vcf] } // channel: [ val(meta), path(vcf) ]
    tbi     = ch_vcf_tbi.map { meta, _vcf, tbi -> [meta, tbi] } // channel: [ val(meta), path(tbi) ]
    depth   = RUN_HIFICNV.out.depth // channel: [ val(meta), path(bw) ]
    copynum = RUN_HIFICNV.out.copynum // channel: [ val(meta), path(bedgraph) ]
    maf     = RUN_HIFICNV.out.maf // channel: [ val(meta), path(bw) ]
}
