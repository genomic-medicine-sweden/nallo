include { MODKIT_PILEUP            } from '../../../modules/nf-core/modkit/pileup/main'
include { TABIX_BGZIPTABIX         } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { MODKIT_BEDMETHYLTOBIGWIG } from '../../../modules/nf-core/modkit/bedmethyltobigwig/main'
include { BEDTOOLS_SORT            } from '../../../modules/nf-core/bedtools/sort/main'

workflow METHYLATION {

    take:
    ch_bam_bai             // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_fai                 // channel: [ val(meta), fai ]
    ch_bed                 // channel: [ val(meta), bed ]
    modcodes               // String or List

    main:
    ch_versions = Channel.empty()

    // Performs pileups per haplotype if the phasing workflow is on, set in config
    MODKIT_PILEUP (ch_bam_bai, ch_fasta, ch_bed)
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)

    BEDTOOLS_SORT (
        MODKIT_PILEUP.out.bed.transpose(),
        []
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    TABIX_BGZIPTABIX ( BEDTOOLS_SORT.out.sorted )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    // Only convert files with content
    BEDTOOLS_SORT.out.sorted
        .filter { _meta, bed -> bed.size() > 0 }
        .set { ch_bedmethyl_to_bigwig_in }

    MODKIT_BEDMETHYLTOBIGWIG ( ch_bedmethyl_to_bigwig_in, ch_fai, modcodes)
    ch_versions = ch_versions.mix(MODKIT_BEDMETHYLTOBIGWIG.out.versions)

    emit:
    bed      = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, bed, _tbi -> [ meta, bed ] } // channel: [ val(meta), path(bed) ]
    tbi      = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, _bed, tbi -> [ meta, tbi ] } // channel: [ val(meta), path(tbi) ]
    bigwig   = MODKIT_BEDMETHYLTOBIGWIG.out.bw
    versions = ch_versions // channel: [ versions.yml ]
}
