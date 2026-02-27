//
// Workflow to call SNVs
//

include { BCFTOOLS_PLUGINFIXPLOIDY                          } from '../../../modules/nf-core/bcftools/pluginfixploidy/main'
include { BEDTOOLS_INTERSECT                                } from '../../../modules/nf-core/bedtools/intersect/main'
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
    variant_caller                  // string: which variant caller to use, e.g. "deepvariant"
    sentieon_tech                   // string: which sequencing tech produced the reads (sentieon)


    main:
    ch_versions   = channel.empty()

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
                meta, _bed ->
                male:   meta.sex == 1
                female: meta.sex == 2
            }
            .set { ch_bed }

        ch_male_diploid_intersect_in   = makeIntersectChannel(ch_sentieon_male_diploid_bed, ch_bed.male, "diploid")
        ch_female_diploid_intersect_in = makeIntersectChannel(ch_sentieon_female_diploid_bed, ch_bed.female, "diploid")
        ch_male_haploid_intersect_in   = makeIntersectChannel(ch_sentieon_male_haploid_bed, ch_bed.male, "haploid")

        ch_male_diploid_intersect_in
            .mix(ch_female_diploid_intersect_in)
            .mix(ch_male_haploid_intersect_in)
            .set { ch_bedtools_intersect_in }

        BEDTOOLS_INTERSECT(
            ch_bedtools_intersect_in,
            [[], []],
        )
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)

        BEDTOOLS_INTERSECT.out.intersect
            .branch {
                meta, intersected_bed ->
                diploid: meta.ploidy == "diploid"
                    [ meta - meta.subMap('ploidy'), intersected_bed  ]
                haploid: meta.ploidy == "haploid"
                    [ meta - meta.subMap('ploidy'), intersected_bed  ]
            }
            .set { ch_intersected_calling_intervals }

        ch_intersected_calling_intervals.haploid
            .map {
                meta, bed ->
                if(bed && bed.size() > 0) {
                    [ meta, bed ]
                } else {
                    [ meta, [] ]
                }
            }
            .mix(
                ch_bed.female.map { meta, _bed -> [ meta, [] ] }
            )
            .set { ch_haploid_regions_out }

        ch_intersected_calling_intervals.diploid
            .set { ch_diploid_regions_out }

        ch_bam_bai
            .join(ch_diploid_regions_out)
            .join(ch_haploid_regions_out)
            .set {
                ch_dnascope_in
            }

        DNASCOPE_LONGREAD(
            ch_dnascope_in,
            ch_fasta,
            ch_fai,
            ch_sentieon_model_bundle,
            sentieon_tech,
        )

        BCFTOOLS_PLUGINFIXPLOIDY(
            DNASCOPE_LONGREAD.out.vcf.join(DNASCOPE_LONGREAD.out.vcf_tbi),
            [],
            [],
            [],
            []
        )

        ch_vcf        = BCFTOOLS_PLUGINFIXPLOIDY.out.vcf
        ch_index      = BCFTOOLS_PLUGINFIXPLOIDY.out.tbi
        ch_gvcf       = DNASCOPE_LONGREAD.out.gvcf
        ch_gvcf_index = DNASCOPE_LONGREAD.out.gvcf_tbi
    }

    emit:
    vcf        = ch_vcf        // channel: [ val(meta), path(vcf) ]
    index      = ch_index      // channel: [ val(meta), path(tbi) ]
    gvcf       = ch_gvcf       // channel: [ val(meta), path(gvcf) ]
    gvcf_index = ch_gvcf_index // channel: [ val(meta), path(tbi) ]
    versions   = ch_versions   // channel: [ path(versions.yml) ]
}

def makeIntersectChannel(ch_sentieon_bed, ch_bed, ploidy_label) {
    ch_sentieon_bed
        .map { _meta, sentieon_regions -> sentieon_regions }
        .combine(ch_bed)
        .map { sentieon_regions, meta, sample_call_regions ->
            [ meta + [ ploidy: ploidy_label ], sample_call_regions, sentieon_regions ]
        }
}
