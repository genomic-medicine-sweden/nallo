include { BCFTOOLS_MERGE                 } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_QUERY                 } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER              } from '../../../modules/nf-core/bcftools/reheader/main'
include { CREATE_SAMPLES_HAPLOTYPES_FILE } from '../../../modules/local/create_samples_haplotypes_file/main'
include { MERGE_JSON                     } from '../../../modules/local/merge_json/main'
include { PARAPHASE                      } from '../../../modules/nf-core/paraphase/main'
include { SAMTOOLS_CONVERT               } from '../../../modules/nf-core/samtools/convert/main'

workflow CALL_PARALOGS {

    take:
    bam_bai     // channel: [ val(meta), bam, bai ]
    fasta       // channel: [ val(meta), fasta    ]
    fai         // channel: [ val(meta), fai      ]
    cram_output // bool: Publish alignments as CRAM (true) or BAM (false)

    main:
    ch_versions = Channel.empty()

    PARAPHASE (
        bam_bai,
        fasta,
        [[],[]]
    )
    ch_versions = ch_versions.mix(PARAPHASE.out.versions)
    PARAPHASE.out.json
        .map { meta, json -> [ [ 'id': meta.family_id ], json ] }
        .groupTuple()
        .set { ch_merge_json_input }

    MERGE_JSON (
        ch_merge_json_input
    )
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

    PARAPHASE.out.vcf
        .transpose()
        .map { meta, vcf ->
            [ [ 'id' : vcf.simpleName, 'family_id': meta.family_id ], vcf, [] ]
        }
        .set { paraphase_vcf_tbis }

    // Get the "sample" name (which is actually the paraphase region, e.g. hba_hba2_hap1) from the VCF
    BCFTOOLS_QUERY (
        paraphase_vcf_tbis,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions)

    // Create rename file for bcftools reheader, e.g. hba_hba2hap1 -> ${sample}_hba_hba2hap1
    CREATE_SAMPLES_HAPLOTYPES_FILE (
        BCFTOOLS_QUERY.out.output
    )
    ch_versions = ch_versions.mix(CREATE_SAMPLES_HAPLOTYPES_FILE.out.versions)

    paraphase_vcf_tbis
        .join( CREATE_SAMPLES_HAPLOTYPES_FILE.out.samples, failOnMismatch:true, failOnDuplicate:true )
        .set { ch_bcftools_reheader_in }

    // Give meta.id as sample name in the VCF
    BCFTOOLS_REHEADER ( ch_bcftools_reheader_in, [[],[]] )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    BCFTOOLS_REHEADER.out.vcf
        .join( BCFTOOLS_REHEADER.out.index, failOnMismatch:true, failOnDuplicate:true )
        .map { meta, vcf, tbi -> [ [ 'id': meta.family_id ], vcf, tbi ] }
        .groupTuple()
        .set { ch_bcftools_merge_in }

    BCFTOOLS_MERGE (
        ch_bcftools_merge_in,
        fasta,
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    emit:
    bam      = PARAPHASE.out.bam                                         // channel: [ val(meta), path(bam) ]
    bai      = PARAPHASE.out.bai                                         // channel: [ val(meta), path(bai) ]
    cram     = cram_output ? SAMTOOLS_CONVERT.out.cram : Channel.empty() // channel: [ val(meta), path(cram) ]
    crai     = cram_output ? SAMTOOLS_CONVERT.out.crai : Channel.empty() // channel: [ val(meta), path(crai) ]
    json     = MERGE_JSON.out.json                                       // channel: [ val(meta), path(json) ]
    vcf      = BCFTOOLS_MERGE.out.vcf                                    // channel: [ val(meta), path(vcfs) ]
    tbi      = BCFTOOLS_MERGE.out.index                                  // channel: [ val(meta), path(tbis) ]
    versions = ch_versions                                               // channel: [ versions.yml ]
}
