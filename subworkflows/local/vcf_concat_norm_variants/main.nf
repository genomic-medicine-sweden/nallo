include { BCFTOOLS_CONCAT                             } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SINGLESAMPLE } from '../../../modules/nf-core/bcftools/norm/main'
include { VCFEXPRESS                                  } from '../../../modules/nf-core/vcfexpress/main'
include { TABIX_BGZIP                                } from '../../../modules/nf-core/tabix/bgzip/main'

//
// Workflow to concatenate and normalize variants
//
workflow VCF_CONCAT_NORM_VARIANTS {
    take:
    ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_fasta                // channel: [mandatory] [ val(meta), path(fasta) ]
    variant_caller          // string: variant caller to tag the variants with, e.g. "deepvariant"
    ch_vcfexpress_prelude   // path: [mandatory] lua file

    main:

    ch_versions = channel.empty()

    BCFTOOLS_CONCAT(
        ch_vcfs.map { meta, vcfs -> [ meta, vcfs, [] ] },
    )

    // Add caller information to meta so vcfexpress can add the FOUND_IN tag based on sv_caller
    BCFTOOLS_CONCAT.out.vcf
        .map { meta, vcf -> [ meta + [ sv_caller: variant_caller ] , vcf ]
        }
        .set { ch_vcfexpress_input }

    VCFEXPRESS (
        ch_vcfexpress_input,
        ch_vcfexpress_prelude
    )

    TABIX_BGZIP {
        VCFEXPRESS.out.vcf
    }

    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)

    // Remove added caller information in meta
    TABIX_BGZIP.out.output
        .map { meta, vcf -> [ meta - meta.subMap('sv_caller'), vcf, [] ]
        }
        .set { ch_bcftools_norm_input }

    BCFTOOLS_NORM_SINGLESAMPLE(
        ch_bcftools_norm_input,
        ch_fasta,
    )

    emit:
    vcf                 = BCFTOOLS_NORM_SINGLESAMPLE.out.vcf                                         // channel: [ val(meta), path(vcf) ]
    index               = BCFTOOLS_NORM_SINGLESAMPLE.out.tbi.mix(BCFTOOLS_NORM_SINGLESAMPLE.out.csi) // channel: [ val(meta), path(tbi/csi) ]
    bcftools_concat_vcf = BCFTOOLS_CONCAT.out.vcf                                                    // channel: [ val(meta), path(vcf) ]
    versions            = ch_versions                                                              // channel: [ path(versions.yml) ]
}
