include { MODKIT_PILEUP            } from '../../../modules/nf-core/modkit/pileup/main'
include { MODKIT_BEDMETHYLTOBIGWIG } from '../../../modules/nf-core/modkit/bedmethyltobigwig/main'

workflow CALL_METHYLATION_MODKIT {
    take:
    ch_bam_bai // channel: [ val(meta), bam, bai ]
    ch_fasta   // channel: [ val(meta), fasta ]
    ch_fai     // channel: [ val(meta), fai ]
    ch_bed     // channel: [ val(meta), bed ]
    modcodes   // String or List

    main:
    ch_versions = channel.empty()

    // Performs pileups per haplotype if the phasing workflow is on, set in config
    MODKIT_PILEUP(
        ch_bam_bai,
        ch_fasta.combine(ch_fai.map { _meta, fai -> fai }),
        ch_bed,
    )
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)

    // Only convert files with content
    MODKIT_PILEUP.out.bedgz
        .transpose()
        .filter { _meta, bed -> bed.size() > 0 }
        .set { ch_bedmethyl_to_bigwig_in }

    MODKIT_BEDMETHYLTOBIGWIG(ch_bedmethyl_to_bigwig_in, ch_fai, modcodes)
    ch_versions = ch_versions.mix(MODKIT_BEDMETHYLTOBIGWIG.out.versions)

    emit:
    bed      = MODKIT_PILEUP.out.bedgz.transpose().map { meta, bed, _tbi -> [meta, bed] } // channel: [ val(meta), path(bed) ]
    tbi      = MODKIT_PILEUP.out.bedgz.transpose().map { meta, _bed, tbi -> [meta, tbi] } // channel: [ val(meta), path(tbi) ]
    bigwig   = MODKIT_BEDMETHYLTOBIGWIG.out.bw
    versions = ch_versions // channel: [ versions.yml ]
}
