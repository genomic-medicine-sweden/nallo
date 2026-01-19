//
// Workflow to call SNVs
//

include { BEDTOOLS_INTERSECT as CREATE_HAPLOID_REGIONS_BED  } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as CREATE_DIPLOID_REGIONS_BED  } from '../../../modules/nf-core/bedtools/intersect/main'
include { DEEPVARIANT_RUNDEEPVARIANT                        } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { DNASCOPE_LONGREAD_CALL_SNVS as DNASCOPE_LONGREAD  } from '../../../modules/local/sentieon/dnascope_longread/main'

workflow CALL_SNVS {
    take:
    ch_bam_bai_bed                  // channel: [mandatory] [ val(meta), path(bam), path(bai), path(call_regions_bed) ]
    ch_fasta                        // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai                          // channel: [mandatory] [ val(meta), path(fai) ]
    ch_par_bed                      // channel: [mandatory] [ val(meta), path(par_bed) ]
    ch_sentieon_model_bundle        // channel: [mandatory] [ val(meta), path(model_bundle) ]
    ch_sentieon_female_diploid_bed  // channel: [mandatory] [ val(meta), path(female_diploid_bed) ]
    ch_sentieon_male_diploid_bed    // channel: [mandatory] [ val(meta), path(male_diploid_bed) ]
    ch_sentieon_male_haploid_bed    // channel: [mandatory] [ val(meta), path(male_haploid_bed) ]
    variant_caller                  // string: which variant caller to use, i.e. "deepvariant"
    sentieon_tech                   // string: which sequencing tech produced the reads (sentieon)


    main:
    ch_versions   = channel.empty()
    ch_vcf        = channel.empty()
    ch_index      = channel.empty()
    ch_gvcf       = channel.empty()
    ch_gvcf_index = channel.empty()

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

    } else if (variant_caller.equals("sentieon")) {

        ch_bam_bai_bed
            .map { meta, bam, bai, _bed ->
                [ meta, bam, bai ]
            }
            .set { ch_bam_bai }

        ch_bam_bai_bed
            .map { meta, _bam, _bai, bed ->
                [ meta, bed ]
            }
            .branch {
                meta, _bam, _bai, _bed ->
                male: meta.sex == 1
                female: meta.sex == 2
            }
            .set { ch_bed }


        ch_bed.male.view()
        ch_sentieon_male_diploid_bed.view()

        ch_sentieon_male_diploid_bed
            .map { _meta, male_diploid_call_regions ->  male_diploid_call_regions }
            .combine(ch_bed.male)
            .set { ch_foo }

        ch_foo.view()

        // ch_foo.map { male_diploid_call_regions, meta, call_regions -> [meta, call_regions, male_diploid_call_regions] }
        //     .set { ch_male_diploid_intersect_in }

    //  ch_male_diploid_intersect_in.view()

        // ch_sentieon_female_diploid_bed
        //     .map { _meta, female_diploid_call_regions -> female_diploid_call_regions }
        //     .combine(ch_bed.female)
        //     .map { female_diploid_call_regions, meta, call_regions -> [meta, call_regions, female_diploid_call_regions] }
        //     .set { ch_female_diploid_intersect_in }

        // ch_male_diploid_intersect_in
        //     .mix(ch_female_diploid_intersect_in)
        //     .set { ch_diploid_intersect_in }

        // CREATE_DIPLOID_REGIONS_BED(
        //     ch_diploid_intersect_in,
        //     [[], []],
        // )
        // ch_versions = ch_versions.mix(CREATE_DIPLOID_REGIONS_BED.out.versions)

        // ch_bam_bai
        //     .join(CREATE_DIPLOID_REGIONS_BED.out.intersect)
        //     .map { meta, bam, bai, _diploid_meta, bed ->
        //         [meta, bam, bai, bed, []]
        //     }
        //     .set {
        //         ch_dnascope_in
        //     }

        // DNASCOPE_LONGREAD(
        //     ch_dnascope_in,
        //     ch_fasta,
        //     ch_fai,
        //     ch_sentieon_model_bundle,
        //     sentieon_tech,
        // )
        // ch_versions = ch_versions.mix(DNASCOPE_LONGREAD.out.versions)

        // ch_vcf        = DNASCOPE_LONGREAD.out.vcf
        // ch_index      = DNASCOPE_LONGREAD.out.vcf_tbi
        // ch_gvcf       = DNASCOPE_LONGREAD.out.gvcf
        // ch_gvcf_index = DNASCOPE_LONGREAD.out.gvcf_tbi
    }


    emit:
    vcf        = ch_vcf        // channel: [ val(meta), path(vcf) ]
    index      = ch_index      // channel: [ val(meta), path(tbi) ]
    gvcf       = ch_gvcf       // channel: [ val(meta), path(gvcf) ]
    gvcf_index = ch_gvcf_index // channel: [ val(meta), path(tbi) ]
    versions   = ch_versions   // channel: [ path(versions.yml) ]
}
