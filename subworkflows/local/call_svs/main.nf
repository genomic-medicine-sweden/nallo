include { ADD_FOUND_IN_TAG    } from '../../../modules/local/add_found_in_tag/main'
include { BCFTOOLS_MERGE      } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_QUERY      } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER   } from '../../../modules/nf-core/bcftools/reheader/main'
include { CREATE_SAMPLES_FILE } from '../../../modules/local/create_samples_file/main'
include { SEVERUS             } from '../../../modules/nf-core/severus/main'
include { SNIFFLES            } from '../../../modules/nf-core/sniffles/main'

workflow CALL_SVS {

    take:
    ch_bam_bai        // channel [mandatory]: [ val(meta), path(bam), path(bai) ]
    ch_fasta          // channel [mandatory]: [ val(meta), path(fasta) ]
    ch_fai            // channel [mandatory]: [ val(meta), path(fai) ]
    sv_caller         //     val [mandatory]: Which caller to use
    ch_tandem_repeats // channel  [optional]: [ val(meta), path(bed) ]
    ch_bed

    main:
    ch_versions     = Channel.empty()

    // Call SVs
    if (sv_caller == "severus") {

        SEVERUS (
            ch_bam_bai.map { meta, bam, bai -> [ meta, bam, bai, [], [], [] ] },
            ch_tandem_repeats
        )
        ch_versions = ch_versions.mix(SEVERUS.out.versions)

        SEVERUS.out.all_vcf
            .set { ch_vcf }

    } else if (sv_caller == "sniffles") {

        SNIFFLES (
            ch_bam_bai,
            ch_fasta,
            ch_tandem_repeats,
            true,
            false
        )
        ch_versions = ch_versions.mix(SNIFFLES.out.versions)

        SNIFFLES.out.vcf
            .set { ch_vcf }
    }

    // Annotate with FOUND_IN tag
    ADD_FOUND_IN_TAG (
        ch_vcf.map { meta, vcf -> [ meta, vcf, [] ] },
        sv_caller
    )
    ch_versions = ch_versions.mix(ADD_FOUND_IN_TAG.out.versions)

    // Get the sample name from the VCF
    // For Sniffles this is hardcoded as SAMPLE and for Severus it's based on the filename
    BCFTOOLS_QUERY (
        ADD_FOUND_IN_TAG.out.vcf.join(ADD_FOUND_IN_TAG.out.csi),
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions)

    // Creates a "vcf_sample_name meta.id" file for bcftools reheader
    CREATE_SAMPLES_FILE ( BCFTOOLS_QUERY.out.output )
    ch_versions = ch_versions.mix(CREATE_SAMPLES_FILE.out.versions)

    ADD_FOUND_IN_TAG.out.vcf
        .join( CREATE_SAMPLES_FILE.out.samples )
        .map { meta, vcf, samples -> [ meta, vcf, [], samples ] }
        .set { ch_bcftools_reheader_in }

    // Give meta.id as sample name in the VCF
    BCFTOOLS_REHEADER ( ch_bcftools_reheader_in, [[],[]] )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    BCFTOOLS_REHEADER.out.vcf
        .join(BCFTOOLS_REHEADER.out.index)
        .map { meta, vcf, tbi -> [ [ 'id': meta.project ], vcf, tbi ] }
        .groupTuple()
        .set { ch_bcftools_merge_in }

    // Merge the files with new sample names
    BCFTOOLS_MERGE ( ch_bcftools_merge_in, ch_fasta, ch_fai, ch_bed )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    emit:
    ch_sv_calls_vcf     = BCFTOOLS_REHEADER.out.vcf   // channel: [ val(meta), path(vcf) ]
    ch_sv_calls_tbi     = BCFTOOLS_REHEADER.out.index // channel: [ val(meta), path(tbi) ]
    ch_multisample_vcf  = BCFTOOLS_MERGE.out.vcf      // channel: [ val(meta), path(vcf) ]
    ch_multisample_tbi  = BCFTOOLS_MERGE.out.index    // channel: [ val(meta), path(tbi) ]
    versions            = ch_versions                 // channel: [ path(versions.yml) ]
}

