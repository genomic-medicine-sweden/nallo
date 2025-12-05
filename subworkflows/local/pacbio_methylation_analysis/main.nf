//
// PacBio Methylation Analysis
//

include { PBCPGTOOLS } from '../../../modules/local/pbcpgtools/main'
include { METHBAT_PROFILE } from '../../../modules/local/methbat/main'
include { METHBAT_COMPARE } from '../../../modules/local/methbat/main'

workflow PACBIO_METHYLATION_ANALYSIS {

    take:
    ch_bam_bai
    ch_reference
    ch_regions

    main:
    ch_versions = Channel.empty()

    //-----------------------------------------------
    // 1️⃣ Run pb-cpg-tools to generate CpG scores
    //-----------------------------------------------
    PBCPGTOOLS(
        ch_bam_bai,
        ch_reference,
        'aligned_bam_to_cpg_scores'
    )
    ch_versions = ch_versions.mix(PBCPGTOOLS.out.versions)

    //-----------------------------------------------
    // 2️⃣ Prepare MethBat input channel
    //-----------------------------------------------
    ch_regions_final = (
        params.regions
            ? Channel.value(file(params.regions))
            : ch_regions.map { it ? file(it) : null }
    )
        .map { r ->
            (r instanceof List && r.isEmpty()) ? null :
            (r == "" ? null : r)
        }
        .ifEmpty { Channel.value(null) }

    ch_methbat_input =
    PBCPGTOOLS.out.bed_and_bw
        .map { meta, bedgz, bw ->
            tuple(meta.id, bedgz, bw)
        }
        .combine(ch_regions_final)

    METHBAT_PROFILE(ch_methbat_input)

    ch_versions = ch_versions.mix(METHBAT_PROFILE.out.versions)

    //-----------------------------------------------
    // 4️⃣ Collect profiles for cohort-level comparison
    //-----------------------------------------------
    ch_cohort_profiles = METHBAT_PROFILE.out.profile
        | collect
        | filter { it.size() > 1 }          // only run compare if ≥2 samples
        | map { profiles ->
            def out = file("cohort_profiles.txt")
            out.withWriter { w ->
                profiles.each { p -> w << p.text }
            }
            return out
        }

    //-----------------------------------------------
    // 5️⃣ Optional MethBat cohort comparison (commented)
    //-----------------------------------------------
    /*
    ch_methbat_compare_out = ch_cohort_profiles
        | map { cohort_file -> tuple(cohort_file, params.baseline, params.compare) }
        | ifEmpty { Channel.empty() }
        | METHBAT_COMPARE

    ch_versions = ch_versions.mix(ch_methbat_compare_out.out.versions)
    */

    emit:
        cpg_scores        = PBCPGTOOLS.out.bed_and_bw
                                .map { meta, bedgz, bw -> tuple(meta, bedgz) }
        methylation_calls = METHBAT_PROFILE.out.profile
        versions          = ch_versions
}

