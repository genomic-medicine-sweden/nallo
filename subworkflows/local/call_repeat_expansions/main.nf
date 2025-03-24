
workflow CALL_REPEAT_EXPANSIONS {

    take:

    main:
    ch_versions         = Channel.empty()
    ch_sample_vcf       = Channel.empty()
    ch_sample_tbi       = Channel.empty()
    ch_family_vcf       = Channel.empty()
    ch_family_tbi       = Channel.empty()
    ch_sample_bam       = Channel.empty()
    ch_sample_bai       = Channel.empty()

if (str_caller == "trgt") {
}

    emit:
    sample_vcf  = ch_sample_vcf // channel: [ val(meta), path(vcf) ]
    sample_tbi  = ch_sample_tbi // channel: [ val(meta), path(tbi) ]
    family_vcf  = ch_family_vcf // channel: [ val(meta), path(vcf) ]
    family_tbi  = ch_family_tbi // channel: [ val(meta), path(tbi) ]
    sample_bam  = ch_sample_bam // channel: [ val(meta), path(bam) ]
    sample_bai  = ch_sample_bai // channel: [ val(meta), path(bai) ]
    versions    = ch_versions   // channel: [ versions.yml ]
}

