include { CAT_FASTQ as CAT_FASTQ_PATERNAL } from '../../../modules/nf-core/cat/fastq/main'
include { CAT_FASTQ as CAT_FASTQ_MATERNAL } from '../../../modules/nf-core/cat/fastq/main'
include { HIFIASM                         } from '../../../modules/nf-core/hifiasm'
include { YAK_COUNT as YAK_PATERNAL       } from '../../../modules/nf-core/yak/count/main'
include { YAK_COUNT as YAK_MATERNAL       } from '../../../modules/nf-core/yak/count/main'
include { GFASTATS                        } from '../../../modules/nf-core/gfastats/main'

workflow ASSEMBLY {
    take:
    ch_reads // channel: [ val(meta), fastq ]

    main:
    ch_versions = Channel.empty()

    if (params.hifiasm_mode == 'hifi-only') {

        ch_reads
            .groupTuple()
            .map { meta, reads -> [meta, reads, []] }
            .set { hifiasm_in }

        HIFIASM(
            hifiasm_in,
            [[], [], []],
            [[], [], []],
            [[], []],
        )
        ch_versions = ch_versions.mix(HIFIASM.out.versions)
    }
    else if (params.hifiasm_mode == 'trio-binning') {
        // Multiple trios with different parents may not work?
        ch_reads
            .groupTuple()
            .map { meta, _reads -> meta }
            .combine(ch_reads.groupTuple().map { meta, reads -> [meta.id, reads] })
            .branch { kid_meta, sample_id, _reads ->
                kid: sample_id == kid_meta.id
                mom: sample_id == kid_meta.maternal_id
                dad: sample_id == kid_meta.paternal_id
            }
            .set { branch_result }

        branch_result.kid.map { meta, _sample_id, reads -> [meta, reads] }.join(branch_result.dad.map { meta, _sample_id, reads -> [meta, reads] }, remainder: true).join(branch_result.mom.map { meta, _sample_id, reads -> [meta, reads] }, remainder: true).branch { kid_meta, kid_reads, dad_reads, mom_reads ->
            is_trio: dad_reads != null && mom_reads != null
            return tuple(kid_meta, kid_reads, dad_reads, mom_reads)
            no_trio: dad_reads == null || mom_reads == null
            return tuple(kid_meta, kid_reads, [], [])
        }.set { ch_samples }

        // We keep using kid_meta for the parental reads, since they need to go together into hifiasm
        trio_kids = ch_samples.is_trio.map { kid_meta, kid_reads, _dad_reads, _mom_reads -> [kid_meta, kid_reads] }
        trio_dads = ch_samples.is_trio.map { kid_meta, _kid_reads, dad_reads, _mom_reads -> [kid_meta, dad_reads] }
        trio_moms = ch_samples.is_trio.map { kid_meta, _kid_reads, _dad_reads, mom_reads -> [kid_meta, mom_reads] }

        // These can be someone elses mom or dad (parents)
        non_trio_kids = ch_samples.no_trio.map { kid_meta, kid_reads, _dad_reads, _mom_reads -> [kid_meta, kid_reads] }
        // Should just be a kid_meta with empty reads
        non_trio_dads = ch_samples.no_trio.map { kid_meta, _kid_reads, dad_reads, _mom_reads -> [kid_meta, dad_reads] }
        non_trio_moms = ch_samples.no_trio.map { kid_meta, _kid_reads, _dad_reads, mom_reads -> [kid_meta, mom_reads] }

        // Parental samples needs to be merged before YAK
        // For now merge every sample, even if there's only one copy of reads
        CAT_FASTQ_PATERNAL(trio_dads)
        ch_versions = ch_versions.mix(CAT_FASTQ_PATERNAL.out.versions)
        CAT_FASTQ_MATERNAL(trio_moms)
        ch_versions = ch_versions.mix(CAT_FASTQ_MATERNAL.out.versions)

        // Then run YAK
        YAK_PATERNAL(CAT_FASTQ_PATERNAL.out.reads)
        ch_versions = ch_versions.mix(YAK_PATERNAL.out.versions)

        YAK_MATERNAL(CAT_FASTQ_MATERNAL.out.reads)
        ch_versions = ch_versions.mix(YAK_MATERNAL.out.versions)

        paternal_yak_or_empty = YAK_PATERNAL.out.yak.concat(non_trio_dads)
        maternal_yak_or_empty = YAK_MATERNAL.out.yak.concat(non_trio_moms)

        trio_kids
            .concat(non_trio_kids)
            .join(paternal_yak_or_empty, failOnMismatch: true, failOnDuplicate: true)
            .join(maternal_yak_or_empty, failOnMismatch: true, failOnDuplicate: true)
            .multiMap { meta, reads, paternal_yak, maternal_yak ->
                reads: [meta, reads, []]
                yak: [meta, paternal_yak, maternal_yak]
            }
            .set { ch_hifiasm_in }

        HIFIASM(
            ch_hifiasm_in.reads,
            ch_hifiasm_in.yak,
            [[], [], []],
            [[], []],
        )
        ch_versions = ch_versions.mix(HIFIASM.out.versions)
    }

    HIFIASM.out.hap1_contigs
        .map { meta, fasta -> [meta + ['haplotype': 1], fasta] }
        .set { ch_gfastats_paternal_in }

    HIFIASM.out.hap2_contigs
        .map { meta, fasta -> [meta + ['haplotype': 2], fasta] }
        .set { ch_gfastats_maternal_in }

    GFASTATS(
        ch_gfastats_paternal_in.mix(ch_gfastats_maternal_in),
        'fasta',
        '',
        '',
        [[], []],
        [[], []],
        [[], []],
        [[], []],
    )
    ch_versions = ch_versions.mix(GFASTATS.out.versions)

    emit:
    assembled_haplotypes = GFASTATS.out.assembly // channel: [ val(meta), path(fasta) ]
    versions             = ch_versions // channel: [ versions.yml ]
}
