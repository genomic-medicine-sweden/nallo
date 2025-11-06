include { BCFTOOLS_PLUGINSPLIT } from '../../../modules/nf-core/bcftools/pluginsplit/main'
include { CREATE_SPLIT_FILE     } from '../../../modules/local/create_split_file/main'

workflow SPLIT_MULTISAMPLE_VCF {
    take:
    ch_vcf       // channel: [ val(meta), path(vcf) ]

    main:
    ch_versions = Channel.empty()

    // Preparing info for splitting
    // Stripping sample IDs from meta for eaiser joining and because it doesn't make sense in single-sample VCFs
    // We also convert sample IDs to list to make transpose work correctly
    ch_vcf
        .map { meta, _vcf -> [ meta - meta.subMap('sample_ids'), meta.sample_ids.toList() ]}
        .transpose()
        .map { meta, sample_id -> [ meta, sample_id, sample_id + '_' + meta.variant_type ]}
        .set { ch_split_info }

    // Convert to format that makes it easy to join and retrieve sample IDs later
    ch_split_info
        .map { meta, sample_id, basename -> [ meta + [basename: basename], sample_id ] }
        .set { ch_split_names }

    ch_split_info
        .collectFile { meta, sample_id, basename -> [
            "${meta.id}_${meta.variant_type}.tsv",
            "${sample_id}\t-\t${basename}"
        ]}
        .map { file ->
            def components = file.simpleName.tokenize('_')
            def meta = [ id: components[0], variant_type: components[1] ]
            return [ meta, file ]
        }
        .set { ch_split_files }

    ch_vcf
        .map { meta, vcf -> [ meta - meta.subMap('sample_ids'), vcf ] }
        .join(ch_split_files, failOnMismatch:true, failOnDuplicate:true)
        .multiMap { meta, vcf, txt ->
            vcf : [ meta, vcf, [] ]
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
        .join(ch_split_names.dump(tag: 'names'), failOnMismatch: true, failOnDuplicate: true)
        .map { meta, file, sample -> [ [ id : sample, family_id : meta.id, variant_type: meta.variant_type], file ] }
        .set { ch_split_vcf }

    emit:
    split_vcf = ch_split_vcf  // channel: [ val(meta), path(vcf) ]
    versions  = ch_versions   // channel: [ path(versions.yml) ]
}
