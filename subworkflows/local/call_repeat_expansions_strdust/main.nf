include { BCFTOOLS_MERGE    } from '../../../modules/nf-core/bcftools/merge/'
include { STRDUST           } from '../../../modules/nf-core/strdust/'
include { HTSLIB_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { VCFEXPRESS        } from '../../../modules/nf-core/vcfexpress/main'

workflow CALL_REPEAT_EXPANSIONS_STRDUST {
    take:
    ch_bam_bai            // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai                // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed                // channel: [mandatory] [ val(meta), path(bed) ]
    ch_vcfexpress_prelude // path: [mandatory] lua file

    main:

    STRDUST(
        ch_bam_bai,
        ch_fasta,
        ch_fai,
        ch_bed,
    )

    VCFEXPRESS(
        STRDUST.out.vcf,
        ch_vcfexpress_prelude,
        'vcf.gz',
    )

    HTSLIB_BGZIPTABIX(
        VCFEXPRESS.out.vcf.map { meta, vcf -> [meta, vcf, [], []] },
        'compress',
        true,
        'vcf',
    )

    ch_bcftools_merge_in = HTSLIB_BGZIPTABIX.out.output
        .join(HTSLIB_BGZIPTABIX.out.index, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi -> [[id: meta.family_id], vcf, tbi] }
        .groupTuple()
        .map { meta, vcfs, tbis -> [meta, vcfs, tbis, []] }

    BCFTOOLS_MERGE(
        ch_bcftools_merge_in,
        ch_fasta.join(ch_fai, failOnMismatch: true, failOnDuplicate: true).collect(),
    )

    emit:
    sample_vcf = VCFEXPRESS.out.vcf          // channel: [ val(meta), path(vcf) ]
    sample_tbi = HTSLIB_BGZIPTABIX.out.index // channel: [ val(meta), path(tbi) ]
    family_vcf = BCFTOOLS_MERGE.out.vcf      // channel: [ val(meta), path(vcf) ]
    family_tbi = BCFTOOLS_MERGE.out.index    // channel: [ val(meta), path(tbi) ]
}
