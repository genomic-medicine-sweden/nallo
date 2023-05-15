// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SNIFFLES_MULTISAMPLE  } from '../../modules/local/sniffles/multisample'
include { SNIFFLES              } from '../../modules/nf-core/sniffles/main'

workflow STRUCTURAL_VARIANT_CALLING {

    take:
    ch_bam_bai // channel: [ val(meta), [[ bam ], [bai]] ]
    ch_snfs
    ch_fasta
    ch_fai
    ch_tandem_repeats

    main:
    
    ch_sv_calls_vcf = Channel.empty()
    ch_versions     = Channel.empty()
    
    meta           = ch_bam_bai.map{ it[0] }
    reference      = ch_fasta.combine(meta).map{[it[2], it[1]]}
    tandem_repeats = ch_tandem_repeats.combine(meta).map{[it[1], it[0]]}
    
    SNIFFLES (ch_bam_bai, reference, tandem_repeats)

    snfs = SNIFFLES.out.snf.map{ it [1] }.concat(ch_snfs.map{ it[1] }).collect().sort{ it.name }.map{ [it] } 
    
    multisample_input = ch_fasta.map{ it[1] }
                        .combine(ch_fai.map{ it[1] })
                        .combine(ch_tandem_repeats)
                        .combine(snfs)
    
    SNIFFLES_MULTISAMPLE( multisample_input )

    ch_versions = ch_versions.mix(SNIFFLES.out.versions)
    ch_versions = ch_versions.mix(SNIFFLES_MULTISAMPLE.out.versions)
    
    ch_sv_calls_vcf = SNIFFLES_MULTISAMPLE.out.vcf

    emit:
    ch_sv_calls_vcf = SNIFFLES.out.vcf             // channel: 
    ch_multisample  = SNIFFLES_MULTISAMPLE.out.vcf //
    versions        = ch_versions                  // channel: [ versions.yml ]
}

