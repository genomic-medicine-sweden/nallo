include { MODKIT_PILEUP                    } from '../../../modules/nf-core/modkit/pileup/main'
include { TABIX_BGZIPTABIX                 } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { MODKIT_BEDMETHYLTOBIGWIG         } from '../../../modules/nf-core/modkit/bedmethyltobigwig/main'
include { PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES } from '../../../modules/nf-core/pbcpgtools/alignedbamtocpgscores/main'
include { METHBAT_PROFILE                  } from '../../../modules/nf-core/methbat/profile/main'

workflow METHYLATION {

    take:
    ch_bam_bai  // channel: [ val(meta), bam, bai ]
    ch_fasta    // channel: [ val(meta), fasta ]
    ch_fai      // channel: [ val(meta), fai ]
    ch_bed      // channel: [ val(meta), bed ]
    modcodes    // String or List
    ch_regions  // channel: [ val(meta), regions ]
    run_methbat // bool: run methbat steps based on preset
    run_modkit  // bool: run modkit steps based on preset

    main:
    ch_versions = channel.empty()

    ch_region_profile        = channel.empty()
    ch_asm_bed               = channel.empty()

    ch_bed                   = channel.empty()
    ch_tbi                   = channel.empty()
    ch_modkit_bigwig         = channel.empty()

    ch_pbcpg_combined_bw     = channel.empty()
    ch_pbcpg_hap1_bw         = channel.empty()
    ch_pbcpg_hap2_bw         = channel.empty()

    if (run_methbat) {
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

        ch_region_profile    = METHBAT_PROFILE.out.region_profile
        ch_asm_bed           = METHBAT_PROFILE.out.asm_bed

        ch_pbcpg_combined_bw = PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.combined_bigwig
        ch_pbcpg_hap1_bw     = PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap1_bigwig
        ch_pbcpg_hap2_bw     = PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES.out.hap2_bigwig
        }

    if (run_modkit) {
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

        ch_bed           = TABIX_BGZIPTABIX.out.gz_tbi.map{ meta, bed, _tbi -> [meta,bed] }
        ch_tbi           = TABIX_BGZIPTABIX.out.gz_tbi.map{ meta, _bed, tbi -> [meta,tbi] }
        ch_modkit_bigwig = MODKIT_BEDMETHYLTOBIGWIG.out.bw
        }

    emit:
    region_profile        = ch_region_profile    // channel: [ val(meta), path(tsv) ]
    asm_bed               = ch_asm_bed           // channel: [ val(meta), path(bed) ]
    bed                   = ch_bed               // channel: [ val(meta), path(bed) ]
    tbi                   = ch_tbi               // channel: [ val(meta), path(tbi) ]
    modkit_bigwig         = ch_modkit_bigwig     // channel: [ val(meta), path(bw) ]
    pbcpg_biwgig_combined = ch_pbcpg_combined_bw // channel: [ val(meta), path(combined.bw) ]
    pbcpg_biwgig_hap1     = ch_pbcpg_hap1_bw     // channel: [ val(meta), path(hap1.bw) ]
    pbcpg_biwgig_hap2     = ch_pbcpg_hap2_bw     // channel: [ val(meta), path(hap2.bw) ]
    versions              = ch_versions          // channel: [ versions.yml ]
}
