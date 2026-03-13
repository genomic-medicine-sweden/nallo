//
// A subworkflow to add the FOUND_IN tag to a vcf (based on caller used)
//

include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'
include { VCFEXPRESS  } from '../../../modules/nf-core/vcfexpress/main'

workflow ADD_FOUND_IN_TAG {
    take:
    ch_vcfexpress_input                   // channel: [mandatory] [ val(meta), path(vcfs) ]
    _ch_variant_caller                     // channel: [mandatory] [ val(meta), path(vcfs) ]
    prelude_content                       // string
    // use  withName: 'VCFEXPRESS' {
    //          ext.args = { -s 'FOUND_IN=return "${meta.sv_caller}"' -e "return true" }
    //      }


    main:

    def lua_file = file("prelude.lua")
    lua_file.text = prelude_content

    //prelude_ch = channel.of(prelude_content)
      //  .collectFile(name: 'prelude.lua')

    VCFEXPRESS (
        ch_vcfexpress_input,
        lua_file
    )

    TABIX_TABIX (
        VCFEXPRESS.out.vcf
    )

    emit:
    vcf      = VCFEXPRESS.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi      = TABIX_TABIX.out.index // channel: [ val(meta), path(tbi) ]
}
