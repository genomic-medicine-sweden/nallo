include { BCFTOOLS_VIEW } from '../../../modules/nf-core/bcftools/view/main'

// Splits mutli-sample VCFs into single-sample VCFs
workflow SPLIT_MULTISAMPLE_VCF {
    take:
    ch_vcf_vartype       // channel: [ val(meta), path(vcf), val(variant_type) ]
    ch_family_to_samples // channel: [ val(meta), val(list_of_sample_ids) ]

    main:
    ch_versions = channel.empty()

    ch_vcf_vartype
        .combine(ch_family_to_samples, by: 0)
        .transpose()
        .map { meta, vcf, variant_type, sample_id ->
            return [[id: sample_id, family_id: meta.id, variant_type: variant_type], vcf, []]
        }
        .set { ch_vcf_prepared }

    BCFTOOLS_VIEW(
        ch_vcf_prepared,
        [],
        [],
        [],
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    BCFTOOLS_VIEW.out.vcf
        .map { meta, vcf ->
            [meta.subMap(['id', 'family_id']), vcf, meta.variant_type]
        }
        .set { ch_split_vcf }

    emit:
    split_vcf = ch_split_vcf // channel: [ val(meta), path(vcf), val(variant_type) ]
    versions = ch_versions   // channel: [ path(versions.yml) ]
}
