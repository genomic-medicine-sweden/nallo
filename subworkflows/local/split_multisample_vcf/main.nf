include { BCFTOOLS_PLUGINSPLIT } from '../../../modules/nf-core/bcftools/pluginsplit/main'

// Splits mutli-sample VCFs into single-sample VCFs
// Requires the meta map to contain a 'sample_ids' key with a set of sample IDs in each VCF
// and a 'variant_type' key with the variant type (e.g. 'snv', 'sv')
workflow SPLIT_MULTISAMPLE_VCF {
    take:
    ch_vcf_vartype       // channel: [ val(meta), path(vcf), val(variant_type) ]
    ch_family_to_samples // channel: [ val(family_id), val(list_of_sample_ids) ]

    main:
    ch_versions = channel.empty()

    // Preparing info for splitting
    // Stripping sample IDs from meta for eaiser joining and because it doesn't make sense in single-sample VCFs
    // We also convert sample IDs to list to make transpose work correctly
    ch_vcf_vartype
        .map { meta, vcf, variant_type -> [ meta.id, meta, vcf, variant_type ] }
        .combine(ch_family_to_samples, by:0) // We can have multiple VCFs per family (e.g. SNVs and SVs)
        .map { _family_id, meta, _vcf, variant_type, sample_ids -> [ meta, variant_type, sample_ids ]}
        .transpose()
        .map { meta, variant_type, sample_id -> [ meta, sample_id, variant_type, sample_id + '_' + variant_type ]}
        .set { ch_split_info }

    // Convert to format that makes it easy to join and retrieve sample IDs later
    ch_split_info
        .map { meta, sample_id, variant_type, basename -> [ meta + [basename: basename, variant_type: variant_type], sample_id ] }
        .set { ch_split_names }

    // Write split info to file to pass to BCFTOOLS_PLUGINSPLIT
    // The dash in the second column indicates that we don't want to rename samples in the output VCFs
    ch_split_info
        .collectFile(newLine: true) { meta, sample_id, variant_type, basename -> [
            "${meta.id}_${variant_type}.tsv",
            "${sample_id}\t-\t${basename}"
        ]}
        .map { file ->
            def components = file.simpleName.tokenize('_')
            def meta = [ id: components[0], variant_type: components[1] ]
            return [ meta, file ]
        }
        .set { ch_split_files }

    // Join VCFs and split files and then split them to ensure correct ordering
    ch_vcf_vartype
        .map { meta, vcf, variant_type -> [ [id : meta.id, variant_type: variant_type ], meta, vcf ] }
        .join(ch_split_files, failOnMismatch:true, failOnDuplicate:true)
        .multiMap { join_key, meta, vcf, txt ->
            vcf : [ meta + [variant_type: join_key.variant_type], vcf, [] ]
            txt : txt
        }
        .set { ch_vcf_prepared }

    BCFTOOLS_PLUGINSPLIT (
        ch_vcf_prepared.vcf,
        ch_vcf_prepared.txt,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_PLUGINSPLIT.out.versions)

    BCFTOOLS_PLUGINSPLIT.out.vcf
        .transpose()
        .map { meta, file -> [ meta + [ basename: file.simpleName ], file ]}
        .join(ch_split_names, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, file, sample -> [ [ id : sample, family_id : meta.id ], file, meta.variant_type ] }
        .set { ch_split_vcf }

    emit:
    split_vcf = ch_split_vcf  // channel: [ val(meta), path(vcf), val(variant_type) ]
    versions  = ch_versions   // channel: [ path(versions.yml) ]
}
