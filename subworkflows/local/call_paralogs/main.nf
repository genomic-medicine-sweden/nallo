include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_PER_FAMILY } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_PER_SAMPLE } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_QUERY                              } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER                           } from '../../../modules/nf-core/bcftools/reheader/main'
include { CREATE_SAMPLES_HAPLOTYPES_FILE              } from '../../../modules/local/create_samples_haplotypes_file/main'
include { MERGE_JSON                                  } from '../../../modules/local/merge_json/main'
include { PARAPHASE                                   } from '../../../modules/nf-core/paraphase/main'
include { SAMTOOLS_CONVERT                            } from '../../../modules/nf-core/samtools/convert/main'

workflow CALL_PARALOGS {

    take:
    bam_bai     // channel: [ val(meta), bam, bai ]
    fasta       // channel: [ val(meta), fasta    ]
    fai         // channel: [ val(meta), fai      ]
    cram_output // bool: Publish alignments as CRAM (true) or BAM (false)

    main:
    ch_versions = Channel.empty()

    PARAPHASE ( bam_bai, fasta, [[],[]] )
    ch_versions = ch_versions.mix(PARAPHASE.out.versions)

    PARAPHASE.out.vcf
        .join( PARAPHASE.out.vcf_index, failOnMismatch:true, failOnDuplicate:true )
        .set { paraphase_vcf_tbis }

    MERGE_JSON ( PARAPHASE.out.json )
    ch_versions = ch_versions.mix(MERGE_JSON.out.versions)

    // Publish bam output as CRAM if requested
    if (cram_output) {
        SAMTOOLS_CONVERT (
            PARAPHASE.out.bam.join(PARAPHASE.out.bai, failOnDuplicate: true, failOnMismatch: true),
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)
    }

    // Get the sample name (GENE_hapX) from the VCF
    BCFTOOLS_QUERY (
        paraphase_vcf_tbis,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions)

    // Creates a "vcf_sample_name meta.id" file for bcftools reheader
    CREATE_SAMPLES_HAPLOTYPES_FILE ( BCFTOOLS_QUERY.out.output )
    ch_versions = ch_versions.mix(CREATE_SAMPLES_HAPLOTYPES_FILE.out.versions)

    PARAPHASE.out.vcf
        .join( CREATE_SAMPLES_HAPLOTYPES_FILE.out.samples, failOnMismatch:true, failOnDuplicate:true )
        .map { meta, vcf, samples -> [ meta, vcf, [], samples ] }
        .set { ch_bcftools_reheader_in }

    // Give meta.id as sample name in the VCF
    BCFTOOLS_REHEADER ( ch_bcftools_reheader_in, [[],[]] )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    BCFTOOLS_REHEADER.out.vcf
        .join( BCFTOOLS_REHEADER.out.index, failOnMismatch:true, failOnDuplicate:true )
        .groupTuple()
        .set { ch_bcftools_merge_in }

    BCFTOOLS_MERGE_PER_SAMPLE(
        ch_bcftools_merge_in,
        fasta,
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE_PER_SAMPLE.out.versions)

    BCFTOOLS_MERGE_PER_SAMPLE.out.vcf
        .join( BCFTOOLS_MERGE_PER_SAMPLE.out.index, failOnMismatch:true, failOnDuplicate:true )
        .map { meta, vcf, tbi -> [ [ 'id': meta.family_id ], vcf, tbi ] }
        .groupTuple()
        .set { bcftools_merge_family_in }

    BCFTOOLS_MERGE_PER_FAMILY ( bcftools_merge_family_in, fasta, [[],[]], [[],[]] )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE_PER_FAMILY.out.versions)




    emit:
    bam      = PARAPHASE.out.bam                     // channel: [ val(meta), path(bam) ]
    bai      = PARAPHASE.out.bai                     // channel: [ val(meta), path(bai) ]
    json     = MERGE_JSON.out.json                   // channel: [ val(meta), path(json) ]
    vcf      = BCFTOOLS_MERGE_PER_FAMILY.out.vcf     // channel: [ val(meta), path(vcfs) ]
    tbi      = BCFTOOLS_MERGE_PER_FAMILY.out.index   // channel: [ val(meta), path(tbis) ]
    versions = ch_versions                           // channel: [ versions.yml ]
}

