include { TRGT_GENOTYPE    } from '../../../modules/nf-core/trgt/genotype/main'
include { SAMTOOLS_INDEX   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT    } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'
include { BCFTOOLS_SORT    } from '../../../modules/nf-core/bcftools/sort/main'
include { TRGT_MERGE       } from '../../../modules/nf-core/trgt/merge/main'
include { VCFEXPRESS       } from '../../../modules/nf-core/vcfexpress/main'

workflow CALL_REPEAT_EXPANSIONS_TRGT {
    take:
    ch_bam_bai              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta                // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai                  // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed                  // channel: [mandatory] [ val(meta), path(bed) ]
    cram_output             // bool: Publish alignments as CRAM (true) or BAM (false)
    ch_vcfexpress_prelude   // channel: [mandatory] [ val(meta), path(lua) ]

    main:

    ch_versions = channel.empty()

    ch_bam_bai
        .map { meta, bam, bai -> [meta, bam, bai, meta.sex == 1 ? 'XY' : 'XX'] }
        .set { ch_trgt_input }

    // Run TRGT
    TRGT_GENOTYPE(
        ch_trgt_input,
        ch_fasta,
        ch_fai,
        ch_bed,
    )

    // Sort and index bam
    SAMTOOLS_SORT(
        TRGT_GENOTYPE.out.bam,
        [[], []],
        '',
    )

    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Publish spanning reads as CRAM if requested
    if (cram_output) {
        SAMTOOLS_CONVERT(
            SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, failOnDuplicate: true, failOnMismatch: true),
            ch_fasta,
            ch_fai,
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)
    }

    _variant_caller = "TRGT"

    TRGT_GENOTYPE.out.vcf
        .multiMap { meta, vcf ->
            vcf: [ meta, vcf ]
            sv_caller: meta.variant_caller
        }
        .set { ch_vcfexpress_input }

    ch_lua_file = ch_vcfexpress_prelude.map { meta, lua -> lua }

    // Add FOUND_IN=TRGT tag
    VCFEXPRESS (
        ch_vcfexpress_input.vcf,
        ch_lua_file
    )

    // Sort and index bcf
    BCFTOOLS_SORT(
        VCFEXPRESS.out.vcf
    )

    // Add sample IDs for all XY samples in family to meta for later repeat annotation with strdrop
    BCFTOOLS_SORT.out.vcf
        .join(BCFTOOLS_SORT.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi -> [[id: meta.family_id], meta, vcf, tbi] }
        .groupTuple()
        .map { meta, sample_metas, vcf, tbi ->
            def xy_ids = sample_metas
                .findAll { sample_meta -> sample_meta.sex == 1 }
                .collect { sample_meta -> sample_meta.id }
                .sort()

            [meta + [xy_samples: xy_ids], vcf, tbi]
        }
        .set { ch_trgt_merge_in }

    TRGT_MERGE(
        ch_trgt_merge_in,
        [[], []],
        [[], []],
    )

    emit:
    sample_vcf = BCFTOOLS_SORT.out.vcf  // channel: [ val(meta), path(vcf) ]
    sample_tbi = BCFTOOLS_SORT.out.tbi  // channel: [ val(meta), path(tbi) ]
    family_vcf = TRGT_MERGE.out.vcf     // channel: [ val(meta), path(vcf) ]
    family_tbi = TRGT_MERGE.out.index   // channel: [ val(meta), path(tbi) ]
    sample_bam = SAMTOOLS_SORT.out.bam  // channel: [ val(meta), path(bam) ]
    sample_bai = SAMTOOLS_INDEX.out.bai // channel: [ val(meta), path(bai) ]
    versions   = ch_versions            // channel: [ versions.yml ]
}
