include { BCFTOOLS_CONCAT                             } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SINGLESAMPLE } from '../../../modules/nf-core/bcftools/norm/main'
include { VCFEXPRESS                                  } from '../../../modules/nf-core/vcfexpress/main'

//
// Workflow to concatenate and normalize variants
//
workflow VCF_CONCAT_NORM_VARIANTS {
    take:
    ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_fasta                // channel: [mandatory] [ val(meta), path(fasta) ]
    variant_caller          // string: variant caller to tag the variants with, e.g. "deepvariant"
    ch_vcfexpress_prelude   // channel: [mandatory] [ val(meta), path(lua) ]

    main:

    BCFTOOLS_CONCAT(
        ch_vcfs.map { meta, vcfs -> [meta, vcfs, []] },
    )

    BCFTOOLS_CONCAT.out.vcf
        .multiMap { meta, vcf ->
            vcf: [ meta, vcf ]
            sv_caller: meta.variant_caller
        }
        .set { ch_vcfexpress_input }

    ch_lua_file = ch_vcfexpress_prelude.map { meta, lua -> lua }

    VCFEXPRESS (
        ch_vcfexpress_input.vcf,
        ch_lua_file
    )

    BCFTOOLS_NORM_SINGLESAMPLE(
        VCFEXPRESS.out.vcf.map { meta, vcf -> [meta, vcf, []] },
        ch_fasta,
    )

    emit:
    vcf                 = BCFTOOLS_NORM_SINGLESAMPLE.out.vcf                                         // channel: [ val(meta), path(vcf) ]
    index               = BCFTOOLS_NORM_SINGLESAMPLE.out.tbi.mix(BCFTOOLS_NORM_SINGLESAMPLE.out.csi) // channel: [ val(meta), path(tbi/csi) ]
    bcftools_concat_vcf = BCFTOOLS_CONCAT.out.vcf                                                    // channel: [ val(meta), path(vcf) ]
}
