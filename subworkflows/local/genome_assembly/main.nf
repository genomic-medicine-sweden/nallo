include { CAT_FASTQ                   } from '../../../modules/nf-core/cat/fastq/main'
include { HIFIASM as HIFIASM_BINS     } from '../../../modules/nf-core/hifiasm'
include { HIFIASM as HIFIASM_ASSEMBLY } from '../../../modules/nf-core/hifiasm'
include { YAK_COUNT                   } from '../../../modules/nf-core/yak/count/main'
include { GFASTATS                    } from '../../../modules/nf-core/gfastats/main'
include { FIND_CONCATENATE            } from '../../../modules/nf-core/find/concatenate/main'


// This subworkflow assembles and outputs haplotypes from a set of reads (grouped per sample), using hifiasm and gfastats.
// It assumes that while each sample can have multiple files, each sample belongs to one family at most.
workflow GENOME_ASSEMBLY {
    take:
    ch_reads // channel: [ val(meta), fastqs ]
    trio_binning // bool: Should we use trio binning mode where possible?
    concat_assemblies // bool: Should we concatenate haplotypes per sample?

    main:
    if (trio_binning) {
        // First, we need to branch the samples based on their relationship
        ch_branched_samples = ch_reads.branch { meta, _reads ->
            def is_parent = meta.relationship in ['father', 'mother']
            paired_parents: is_parent && meta.has_other_parent
            children_with_both_parents: meta.relationship == 'child' && meta.two_parents
            other: true
        }

        // Then, the files from parents of children with both parents will need to be concatenated before yak
        // in case there are multiple files for the same parent.
        ch_paired_parents_for_yak = ch_branched_samples.paired_parents.branch { _meta, fastqs ->
            cat: fastqs.size() > 1
            no_cat: fastqs.size() == 1
        }

        CAT_FASTQ(
            ch_paired_parents_for_yak.cat
        )

        YAK_COUNT(
            CAT_FASTQ.out.reads.concat(ch_paired_parents_for_yak.no_cat)
        )

        ch_yak_output = YAK_COUNT.out.yak
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

        // Creates the input for trio-binned assemblies (children with both parents)
        ch_with_both_parents = ch_branched_samples.children_with_both_parents
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(ch_yak_output.paternal)
            .join(ch_yak_output.maternal)
            .map { _id, meta, reads, yak_paternal, yak_maternal ->
                [meta, reads, yak_paternal, yak_maternal]
            }

        // Create the input for hifiasm by combining the non-trio binned samples with the trio-binned samples.
        ch_hifiasm_in = ch_branched_samples.other
            .concat(ch_branched_samples.paired_parents)
            .map { meta, fastqs ->
                [meta, fastqs, [], []]
            }
            .concat(ch_with_both_parents)
            .multiMap { meta, reads, yak_paternal, yak_maternal ->
                reads: [meta, reads, []]
                yak: [meta, yak_paternal, yak_maternal]
            }
    }
    else {
        ch_hifiasm_in = ch_reads.multiMap { meta, reads ->
            reads: [meta, reads, []]
            yak: [meta, [], []]
        }
    }

    HIFIASM_BINS(
        ch_hifiasm_in.reads,
        ch_hifiasm_in.yak,
        [[], [], []],
        [[], []],
    )

    // Explicitly key bins/reads/yak by sample ID before assembly so each sample gets its own bins and yaks.
    ch_hifiasm_assembly_in = ch_hifiasm_in.reads
        .join(ch_hifiasm_in.yak, failOnMismatch: true, failOnDuplicate: true)
        .join(HIFIASM_BINS.out.bin_files, failOnMismatch: true, failOnDuplicate: true)
        .multiMap { meta, reads, ul_reads, yak_paternal, yak_maternal, bin_files ->
            reads: [meta, reads, ul_reads]
            bins: [meta, bin_files]
            yak: [meta, yak_paternal, yak_maternal]
        }

    HIFIASM_ASSEMBLY(
        ch_hifiasm_assembly_in.reads,
        ch_hifiasm_assembly_in.yak,
        [[], [], []],
        ch_hifiasm_assembly_in.bins,
    )

    ch_gfastats_paternal_in = HIFIASM_ASSEMBLY.out.hap1_contigs.map { meta, fasta -> [meta + ['haplotype': 1], fasta] }

    ch_gfastats_maternal_in = HIFIASM_ASSEMBLY.out.hap2_contigs.map { meta, fasta -> [meta + ['haplotype': 2], fasta] }

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

    if (concat_assemblies) {
        ch_assemblies_to_concatenate = GFASTATS.out.assembly
            .map { meta, bam -> [meta - meta.subMap('haplotype'), bam] }
            .groupTuple(size: 2)

        FIND_CONCATENATE(
            ch_assemblies_to_concatenate
        )

        ch_concatenated_haplotypes = FIND_CONCATENATE.out.file_out
    }
    else {
        ch_concatenated_haplotypes = channel.empty()
    }

    emit:
    assembled_haplotypes    = GFASTATS.out.assembly // channel: [ val(meta), path(fasta) ]
    assembly_summary        = GFASTATS.out.assembly_summary // channel: [ val(meta), path(assembly_summary) ]
    concatenated_haplotypes = ch_concatenated_haplotypes // channel: [ val(meta), path(fasta) ]
}
