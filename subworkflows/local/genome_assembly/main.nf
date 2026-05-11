include { CAT_FASTQ                   } from '../../../modules/nf-core/cat/fastq/main'
include { HIFIASM as HIFIASM_BINS     } from '../../../modules/nf-core/hifiasm'
include { HIFIASM as HIFIASM_ASSEMBLY } from '../../../modules/nf-core/hifiasm'
include { YAK_COUNT                   } from '../../../modules/nf-core/yak/count/main'
include { GFASTATS                    } from '../../../modules/nf-core/gfastats/main'

// This subworkflow assembles and outputs haplotypes from a set of reads (grouped per sample), using hifiasm and gfastats.
// It assumes that while each sample can have multiple files, each sample belongs to one family at most.
workflow GENOME_ASSEMBLY {
    take:
    ch_reads // channel: [ val(meta), fastqs ]
    trio_binning // bool: Should we use trio binning mode where possible?

    main:
    if (trio_binning) {
        // First, we need to branch the samples based on their relationship
        ch_reads
            .branch { meta, _reads ->
                def is_parent = meta.relationship in ['father', 'mother']
                paired_parents: is_parent && meta.has_other_parent
                children_with_both_parents: meta.relationship == 'child' && meta.two_parents
                other: true
            }
            .set { ch_branched_samples }

        // Then, the files from parents of children with both parents will need to be concatenated before yak
        // in case there are multiple files for the same parent.
        ch_branched_samples.paired_parents
            .branch { _meta, fastqs ->
                cat: fastqs.size() > 1
                no_cat: fastqs.size() == 1
            }
            .set { ch_paired_parents_for_yak }

        CAT_FASTQ(
            ch_paired_parents_for_yak.cat
        )

        YAK_COUNT(
            CAT_FASTQ.out.reads.concat(ch_paired_parents_for_yak.no_cat)
        )

        YAK_COUNT.out.yak
            .flatMap { meta, yak ->
                (meta.children ?: []).collect { child_id ->
                    [child_id, meta, yak]
                }
            }
            .branch { child_id, meta, yak ->
                paternal: meta.relationship == 'father'
                return [child_id, yak]
                maternal: meta.relationship == 'mother'
                return [child_id, yak]
            }
            .set { ch_yak_output }

        // Creates the input for trio-binned assemblies (children with both parents)
        ch_branched_samples.children_with_both_parents
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(ch_yak_output.paternal)
            .join(ch_yak_output.maternal)
            .map { _id, meta, reads, yak_paternal, yak_maternal ->
                [meta, reads, yak_paternal, yak_maternal]
            }
            .set { ch_with_both_parents }

        // Create the input for hifiasm by combining the non-trio binned samples with the trio-binned samples.
        ch_branched_samples.other
            .concat(ch_branched_samples.paired_parents)
            .map { meta, fastqs ->
                [meta, fastqs, [], []]
            }
            .concat(ch_with_both_parents)
            .multiMap { meta, reads, yak_paternal, yak_maternal ->
                reads: [meta, reads, []]
                yak: [meta, yak_paternal, yak_maternal]
            }
            .set { ch_hifiasm_in }
    }
    else {
        ch_reads
            .multiMap { meta, reads ->
                reads: [meta, reads, []]
                yak: [meta, [], []]
            }
            .set { ch_hifiasm_in }
    }

    HIFIASM_BINS(
        ch_hifiasm_in.reads,
        ch_hifiasm_in.yak,
        [[], [], []],
        [[], []],
    )

    // Explicitly key bins/reads/yak by sample ID before assembly so each sample gets its own bins and yaks.
    ch_hifiasm_in.reads
        .join(ch_hifiasm_in.yak, failOnMismatch: true, failOnDuplicate: true)
        .join(HIFIASM_BINS.out.bin_files, failOnMismatch: true, failOnDuplicate: true)
        .multiMap { meta, reads, ul_reads, yak_paternal, yak_maternal, bin_files ->
            reads: [meta, reads, ul_reads]
            bins: [meta, bin_files]
            yak: [meta, yak_paternal, yak_maternal]
        }
        .set { ch_hifiasm_assembly_in }

    HIFIASM_ASSEMBLY(
        ch_hifiasm_assembly_in.reads,
        ch_hifiasm_assembly_in.yak,
        [[], [], []],
        ch_hifiasm_assembly_in.bins,
    )

    HIFIASM_ASSEMBLY.out.hap1_contigs
        .map { meta, fasta -> [meta + ['haplotype': 1], fasta] }
        .set { ch_gfastats_paternal_in }

    HIFIASM_ASSEMBLY.out.hap2_contigs
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

    emit:
    assembled_haplotypes = GFASTATS.out.assembly // channel: [ val(meta), path(fasta) ]
    assembly_summary = GFASTATS.out.assembly_summary // channel: [ val(meta), path(assembly_summary) ]
}
