include { BCFTOOLS_PLUGINSPLIT } from '../../../modules/nf-core/bcftools/pluginsplit/main'
include { CREATE_SPLIT_FILE     } from '../../../modules/local/create_split_file/main'

workflow SPLIT_MULTISAMPLE_VCF {
    take:
    ch_vcf       // channel: [ val(meta), path(vcf) ]
    suffix       // string: suffix for split file naming (e.g., "_snv" or "_sv")

    main:
    ch_versions = Channel.empty()

    // Create sample files for bcftools +split
    ch_vcf
        .map { meta, _vcf -> [ meta, meta.sample_ids ] }
        .set { ch_create_split_file_in }

    CREATE_SPLIT_FILE ( ch_create_split_file_in, suffix )

    // Collect the files' content into channels
    CREATE_SPLIT_FILE.out.txt
        .splitCsv(sep: "\t", elem: 1, header: ['sample', 'dash', 'basename'])
        .map { meta, row -> [ meta + [ basename: row.basename ], row.sample ] }
        .set { ch_split_names } // [ [region, family, samples, basename], sample ]

    ch_vcf
        .join(CREATE_SPLIT_FILE.out.txt, failOnMismatch:true, failOnDuplicate:true)
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
        .join(ch_split_names, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, file, sample -> [ [ id : sample, family_id : meta.id ], file ] }
        .set { ch_split_vcf }

    emit:
    split_vcf = ch_split_vcf  // channel: [ val(meta), path(vcf) ]
    versions  = ch_versions   // channel: [ path(versions.yml) ]
}
