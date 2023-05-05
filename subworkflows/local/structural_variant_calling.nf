// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SNIFFLES_SINGLESAMPLE      } from '../../modules/local/sniffles/singlesample'
include { SNIFFLES_MULTISAMPLE      } from '../../modules/local/sniffles/multisample'

workflow STRUCTURAL_VARIANT_CALLING {

    take:
    ch_bam_bai // channel: [ val(meta), [[ bam ], [bai]] ]
    ch_snfs
    ch_fasta

    main:
    
    ch_sv_calls_vcf = Channel.empty()
    ch_versions     = Channel.empty()
    
    SNIFFLES_SINGLESAMPLE( ch_bam_bai.combine(ch_fasta.map {it [1] }) )
    SNIFFLES_MULTISAMPLE( SNIFFLES_SINGLESAMPLE.out.sv_snf.map { it [1] }.concat(ch_snfs.map { it[1] }).collect().sort { it.name } )

    ch_versions = ch_versions.mix(SNIFFLES_SINGLESAMPLE.out.versions)
    ch_versions = ch_versions.mix(SNIFFLES_MULTISAMPLE.out.versions)
    
    ch_sv_calls_vcf = SNIFFLES_MULTISAMPLE.out.multisample_vcf


    emit:
    ch_sv_calls_vcf
    
    versions = ch_versions                  // channel: [ versions.yml ]
}

