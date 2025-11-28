include { MODKIT_PILEUP                    } from '../../../modules/nf-core/modkit/pileup/main'
include { TABIX_BGZIPTABIX                 } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { MODKIT_BEDMETHYLTOBIGWIG         } from '../../../modules/nf-core/modkit/bedmethyltobigwig/main'
include { PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES } from '../../../modules/nf-core/pbcpgtools/alignedbamtocpgscores/main'
include { METHBAT_PROFILE                  } from '../../../modules/nf-core/methbat/profile/main'

workflow METHYLATION {

    take:
    ch_bam_bai             // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_fai                 // channel: [ val(meta), fai ]
    ch_bed                 // channel: [ val(meta), bed ]
    modcodes               // String or List
    ch_regions             // channel: [ val(meta), regions ]
    // add boolean based on ONT/Pacbio

    main:
    ch_versions = channel.empty()

    PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES ( ch_bam_bai )
    ch_versions = ch_versions.mix(PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.versions)

    PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.combined_bed
        .mix(PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.combined_bed_index,
             PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap1_bed,
             PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap1_bed_index,
             PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap2_bed,
             PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap2_bed_index)
        .groupTuple()
        .set { ch_methbat_profile_in }

    METHBAT_PROFILE ( ch_methbat_profile_in, ch_regions ) 
    ch_versions = ch_versions.mix(METHBAT_PROFILE.out.versions)

    // Performs pileups per haplotype if the phasing workflow is on, set in config
    MODKIT_PILEUP (ch_bam_bai, ch_fasta, ch_bed)
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)

    MODKIT_PILEUP.out.bed
        .transpose()
        .set { ch_bedmethyl }

    TABIX_BGZIPTABIX ( ch_bedmethyl )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    // Only convert files with content
    ch_bedmethyl
        .filter { _meta, bed -> bed.size() > 0 }
        .set { ch_bedmethyl_to_bigwig_in }

    MODKIT_BEDMETHYLTOBIGWIG ( ch_bedmethyl_to_bigwig_in, ch_fai, modcodes)
    ch_versions = ch_versions.mix(MODKIT_BEDMETHYLTOBIGWIG.out.versions)

    emit:
    region_profile        = METHBAT_PROFILE.out.region_profile                                   // channel: [ val(meta), path(tsv) ]
    asm_bed               = METHBAT_PROFILE.out.asm_bed                                          // channel: [ val(meta), path(bed) ]
    bed                   = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, bed, _tbi -> [ meta, bed ] } // channel: [ val(meta), path(bed) ]
    tbi                   = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, _bed, tbi -> [ meta, tbi ] } // channel: [ val(meta), path(tbi) ]
    modkit_bigwig         = MODKIT_BEDMETHYLTOBIGWIG.out.bw                                      // channel: [ val(meta), path(bw) ]
    pbcpg_biwgig_combined = PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.combined_bigwig                 // channel: [ val(meta), path(combined.bw) ]
    pbcpg_biwgig_hap1     = PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap1_bigwig                     // channel: [ val(meta), path(hap1.bw) ]
    pbcpg_biwgig_hap2     = PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap2_bigwig                     // channel: [ val(meta), path(hap2.bw) ]
    versions              = ch_versions                                                          // channel: [ versions.yml ]
}
