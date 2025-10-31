//
// Workflow to call SNVs
//
include { DEEPVARIANT_RUNDEEPVARIANT } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'

workflow CALL_SNVS {
    take:
    ch_bam_bai_bed // channel: [mandatory] [ val(meta), path(bam), path(bai), path(call_regions_bed) ]
    ch_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai         // channel: [mandatory] [ val(meta), path(fai) ]
    ch_par_bed     // channel: [mandatory] [ val(meta), path(par_bed) ]
    variant_caller //  string: which variant caller to use, i.e. "deepvariant"

    main:
    ch_versions   = Channel.empty()
    ch_vcf        = Channel.empty()
    ch_index      = Channel.empty()
    ch_gvcf       = Channel.empty()
    ch_gvcf_index = Channel.empty()

    if (variant_caller.equals("deepvariant")) {

        DEEPVARIANT_RUNDEEPVARIANT(
            ch_bam_bai_bed,
            ch_fasta,
            ch_fai,
            [[], []],
            ch_par_bed,
        )
        ch_versions = ch_versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions)

        ch_vcf        = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
        ch_index      = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index
        ch_gvcf       = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
        ch_gvcf_index = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index
    }

    emit:
    vcf        = ch_vcf        // channel: [ val(meta), path(vcf) ]
    index      = ch_index      // channel: [ val(meta), path(tbi) ]
    gvcf       = ch_gvcf       // channel: [ val(meta), path(gvcf) ]
    gvcf_index = ch_gvcf_index // channel: [ val(meta), path(tbi) ]
    versions   = ch_versions   // channel: [ path(versions.yml) ]
}
