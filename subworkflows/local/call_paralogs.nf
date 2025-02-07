include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_PER_CASE } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_PER_SAMPLE} from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_QUERY } from '../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER } from '../../modules/nf-core/bcftools/reheader/main'
include { CREATE_SAMPLES_FILE } from '../../modules/local/create_samples_file/main'
include { PARAPHASE } from '../../modules/nf-core/paraphase/main'

workflow CALL_PARALOGS {

    take:
    bam_bai // channel: [ val(meta), bam, bai ]
    fasta   // channel: [ val(meta), fasta ]

    main:
    ch_versions = Channel.empty()

    PARAPHASE ( bam_bai, fasta, [[],[]] )
    ch_versions = ch_versions.mix(PARAPHASE.out.versions)

    PARAPHASE.out.vcf
        .join(PARAPHASE.out.vcf_index)
        .set { paraphase_vcf_tbis }

    // Get the sample name (GENE_hapX)from the VCF
    BCFTOOLS_QUERY (
        paraphase_vcf_tbis,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions)

    // Creates a "vcf_sample_name meta.id" file for bcftools reheader
    CREATE_SAMPLES_FILE ( BCFTOOLS_QUERY.out.output )
    ch_versions = ch_versions.mix(CREATE_SAMPLES_FILE.out.versions)

    PARAPHASE.out.vcf
        .join( CREATE_SAMPLES_FILE.out.samples )
        .map { meta, vcf, samples -> [ meta, vcf, [], samples ] }
        .set { ch_bcftools_reheader_in }

    // Give meta.id as sample name in the VCF
    BCFTOOLS_REHEADER ( ch_bcftools_reheader_in, [[],[]] )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    BCFTOOLS_REHEADER.out.vcf
        .join( BCFTOOLS_REHEADER.out.index )
        .set { ch_bcftools_merge_in }

    BCFTOOLS_MERGE_PER_SAMPLE ( ch_bcftools_merge_in, [], [], [] )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE_PER_SAMPLE.out.versions)

    BCFTOOLS_MERGE_PER_SAMPLE.out.vcf
        .join( BCFTOOLS_MERGE_PER_SAMPLE.out.index )
        .map{ meta, vcf, tbi -> [ [ 'id': meta.family_id ], vcf, tbi ] }
        .set {bcftools_merge_case_in}

    BCFTOOLS_MERGE_PER_CASE ( bcftools_merge_case_in, [], [], [] )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE_PER_CASE.out.versions)




    emit:
    bam      = PARAPHASE.out.bam                     // channel: [ val(meta), path(bam) ]
    bai      = PARAPHASE.out.bai                     // channel: [ val(meta), path(bai) ]
    json     = PARAPHASE.out.json                    // channel: [ val(meta), path(json) ]
    vcf      = BCFTOOLS_MERGE_PER_CASE.out.vcf       // channel: [ val(meta), path(vcfs) ]
    tbi      = BCFTOOLS_MERGE_PER_CASE.out.index     // channel: [ val(meta), path(tbis) ]
    versions = ch_versions                           // channel: [ versions.yml ]
}

