include { BCFTOOLS_MERGE                 } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_QUERY                 } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER              } from '../../../modules/nf-core/bcftools/reheader/main'
include { GAWK                           } from '../../../modules/nf-core/gawk/main'
include { PARAPHASE                      } from '../../../modules/nf-core/paraphase/main'
include { SAMTOOLS_CONVERT               } from '../../../modules/nf-core/samtools/convert/main'
include { findKeysForValue } from '../utils_nfcore_nallo_pipeline/main.nf'
workflow CALL_PARALOGS {
    take:
    bam_bai     // channel: [ val(meta), bam, bai ]
    fasta       // channel: [ val(meta), fasta ]
    fai         // channel: [ val(meta), fai ]
    cram_output // bool: Publish alignments as CRAM (true) or BAM (false)

    main:
    PARAPHASE(
        bam_bai,
        fasta,
        [[], []],
    )

    PARAPHASE.out.vcf
        .transpose()
        .map { meta, vcf ->
            [['id': vcf.simpleName, 'family_id': meta.family_id], vcf, []]
        }
        .set { paraphase_vcf_tbis }

    // Extract the Paraphase locus identifier from the VCF (e.g. hba_hba2hap1). This is encoded in the VCF as the sample name.
    BCFTOOLS_QUERY(
        paraphase_vcf_tbis,
        [],
        [],
        [],
    )

    /*
     * Create a bcftools reheader mapping file to make Paraphase VCF sample names globally unique.
     *
     * We add the biological sample name as a prefix to the paraphase identifier (e.g. hba_hba2hap1 -> ${sample}_hba_hba2hap1), since bcftools merge requires all sample names across input VCFs to be unique.
     */
    GAWK(
        BCFTOOLS_QUERY.out.output,
        [],
        false,
    )

    paraphase_vcf_tbis
        .join(GAWK.out.output, failOnMismatch: true, failOnDuplicate: true)
        .set { ch_bcftools_reheader_in }

    BCFTOOLS_REHEADER(ch_bcftools_reheader_in, [[], []])

    BCFTOOLS_REHEADER.out.vcf
        .join(BCFTOOLS_REHEADER.out.index, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi -> [meta.family_id, vcf, tbi] }
        .groupTuple()
        .map { family_id, vcfs, tbis -> [['id': family_id], vcfs, tbis, []] }
        .set { ch_reheadered_vcf_tbis_per_family }

    BCFTOOLS_MERGE(
        ch_reheadered_vcf_tbis_per_family,
        fasta.join(fai, failOnMismatch: true, failOnDuplicate: true),
    )

    if (cram_output) {
        SAMTOOLS_CONVERT(
            PARAPHASE.out.bam.join(PARAPHASE.out.bai, failOnDuplicate: true, failOnMismatch: true),
            fasta.join(fai, failOnDuplicate: true, failOnMismatch: true).collect(),
        )
    }

    emit:
    bam      = PARAPHASE.out.bam                                         // channel: [ val(meta), path(bam) ]
    bai      = PARAPHASE.out.bai                                         // channel: [ val(meta), path(bai) ]
    cram     = cram_output ? SAMTOOLS_CONVERT.out.cram : channel.empty() // channel: [ val(meta), path(cram) ]
    crai     = cram_output ? SAMTOOLS_CONVERT.out.crai : channel.empty() // channel: [ val(meta), path(crai) ]
    json     = PARAPHASE.out.json                                        // channel: [ val(meta), path(json) ]
    vcf      = BCFTOOLS_MERGE.out.vcf                                    // channel: [ val(meta), path(vcfs) ]
    tbi      = BCFTOOLS_MERGE.out.index                                  // channel: [ val(meta), path(tbis) ]
}
