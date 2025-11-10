include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main'
include { WHATSHAP_HAPLOTAG } from '../../../modules/local/whatshap/haplotag/main'
include { WHATSHAP_PHASE    } from '../../../modules/local/whatshap/phase/main'

workflow WHATSHAP {
    take:
    ch_snv_vcf   // channel: [ val(meta), path(vcf) ]
    ch_bam_bai   // channel: [ val(meta), path(bam), path(bai) ]
    fasta        // channel: [ val(meta), path(fasta) ]
    fai          // channel: [ val(meta), path(fai) ]

    main:
    ch_versions = Channel.empty()

    // Fix metadata to group by family
    ch_bam_bai
        .map { meta, bam, bai -> [ [ id : meta.family_id ], bam, bai ] }
        .groupTuple()
        .set { ch_bam_bai_grouped }

    ch_snv_vcf
        .map { meta, vcf -> [ [ id: meta.id ], vcf ] }
        .join(ch_bam_bai_grouped, failOnMismatch: true, failOnDuplicate: true)
        .set { ch_whatshap_phase_in}

    WHATSHAP_PHASE(
        ch_whatshap_phase_in,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions)

    WHATSHAP_PHASE.out.vcf_tbi
        .multiMap { meta, vcf, tbi ->
            vcf : [ meta, vcf ]
            tbi : [ meta, tbi ]
        }
        .set { ch_snv_outputs }

    // We cannot use the grouped BAM channel here because WhatsHap can haplotag only one BAM at a time.
    // Using combine instead of join because the VCFs are family-level, not sample-level
    // Therefore, there might be multiple BAMs per VCF and join only keeps the first match
    // (unless failOnDuplicate is true, then we get an error)
    ch_bam_bai
        .map { meta, bam, bai -> [ [ id : meta.family_id ], meta, bam, bai ] }
        .combine(WHATSHAP_PHASE.out.vcf_tbi, by: 0)
        .map { _family_meta, sample_meta, bam, bai, vcf, tbi -> [ sample_meta, vcf, tbi, bam, bai ] }
        .set { ch_whatshap_haplotag_in }

    WHATSHAP_HAPLOTAG (
        ch_whatshap_haplotag_in,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)

    SAMTOOLS_INDEX (
        WHATSHAP_HAPLOTAG.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    WHATSHAP_HAPLOTAG.out.bam
        .join( SAMTOOLS_INDEX.out.bai, failOnMismatch:true, failOnDuplicate:true )
        .set { ch_bam_bai_haplotagged }

    emit:
    phased_family_snvs     = ch_snv_outputs.vcf  // channel: [ val(meta), path(vcf) ]
    phased_family_snvs_tbi = ch_snv_outputs.tbi  // channel: [ val(meta), path(tbi) ]
    haplotagged_bam_bai    = ch_bam_bai_haplotagged     // channel: [ val(meta), path(bam), path(bai) ]
    versions               = ch_versions                // channel: [ path(versions.yml) ]
}
