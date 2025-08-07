include { CAT_FASTQ } from '../../../modules/nf-core/cat/fastq/main'
include { HIFIASM   } from '../../../modules/nf-core/hifiasm'
include { YAK_COUNT } from '../../../modules/nf-core/yak/count/main'
include { GFASTATS  } from '../../../modules/nf-core/gfastats/main'

// This subworkflow assembles and outputs haplotypes from a set of reads, using hifiasm and gfastats.
// It assumes that no sample exists twice in different families.
workflow GENOME_ASSEMBLY {

    take:
    ch_reads      // channel: [ val(meta), fastq ]
    trio_binning  //    bool: Should we use trio binning mode where possible?

    main:
    ch_versions = Channel.empty()

    if (!trio_binning) {
        ch_reads
            .map { meta, fastq ->
                [ groupKey(meta, meta.n_files), fastq ]
            }
            .groupTuple()
            .multiMap { meta, reads ->
                reads : [ meta, reads, [] ]
                yak   : [ [], [], [] ]
            }
            .set { ch_hifiasm_in }

    } else {
        // First, we need to branch the samples based on their relationship
        ch_reads
            .branch { meta, _reads ->
                def is_parent = meta.relationship in ['father', 'mother']
                paired_parents             : is_parent && meta.has_other_parent
                unpaired_parents           : is_parent && !meta.has_other_parent
                children_with_both_parents : meta.relationship == 'child' && meta.two_parents
                children_with_single_parent: meta.relationship == 'child' && !meta.two_parents
                unknown                    : meta.relationship == 'unknown'
            }
            .set { ch_branched_samples }

        // Then, the files from parents of children with both parents will need to be concatenated before yak
        // in case there are multiple files for the same parent.
        ch_branched_samples.paired_parents
            .map { meta, fastq ->
                [ groupKey(meta, meta.n_files), fastq ]
            }
            .groupTuple()
            .branch { _meta, fastqs ->
                cat: fastqs.size() > 1
                no_cat: fastqs.size() == 1
            }
            .set { ch_paired_parents_for_yak }

        CAT_FASTQ (
            ch_paired_parents_for_yak.cat
        )
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

        YAK_COUNT (
            CAT_FASTQ.out.reads.concat(ch_paired_parents_for_yak.no_cat)
        )
        ch_versions = ch_versions.mix(YAK_COUNT.out.versions)

        YAK_COUNT.out.yak
            // Because a parent can have multiple children, and meta.children is a list of all children,
            // we need to return one tuple per child.
            .flatMap { meta, yak ->
                (meta.children ?: []).collect { child_id ->
                    [child_id, meta, yak]
                }
            }
            .branch { child_id, meta, yak ->
                paternal: meta.relationship == 'father'
                    return [ child_id, yak ]
                maternal: meta.relationship == 'mother'
                    return [ child_id, yak ]
            }
            .set { ch_yak_output }

        // Creates the input for trio-binned assemblies (children with both parents)
        ch_branched_samples.children_with_both_parents
            .map { meta, fastq ->
                [ groupKey(meta, meta.n_files), fastq ]
            }
            .groupTuple()
            .map { meta, reads -> [ meta.id, meta, reads ] }
            .join(ch_yak_output.paternal)
            .join(ch_yak_output.maternal)
            .map { _id, meta, reads, yak_paternal, yak_maternal ->
                [ meta, reads, yak_paternal, yak_maternal ]
            }
            .set { ch_with_both_parents }

        // Creates the input for the non-trio binned assemblies
        // (children with a single parent, parents, and unknown samples)
        ch_branched_samples.children_with_single_parent
            .concat(ch_branched_samples.paired_parents)
            .concat(ch_branched_samples.unpaired_parents)
            .concat(ch_branched_samples.unknown)
            .map { meta, fastq ->
                [ groupKey(meta, meta.n_files), fastq ]
            }
            .groupTuple()
            .map { meta, fastqs ->
                [ meta, fastqs, [], [] ]
            }
            .set { ch_other_cases }

        // Combines the two, and creates the input for hifiasm
        ch_with_both_parents
            .concat(ch_other_cases)
            .multiMap { meta, reads, yak_paternal, yak_maternal ->
                reads : [ meta, reads, []                  ]
                yak   : [ meta, yak_paternal, yak_maternal ]
            }
            .set { ch_hifiasm_in }
    }

    HIFIASM (
        ch_hifiasm_in.reads,
        ch_hifiasm_in.yak,
        [[],[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(HIFIASM.out.versions)

    HIFIASM.out.hap1_contigs
        .map { meta, fasta -> [ meta + [ 'haplotype': 1 ], fasta ] }
        .set { ch_gfastats_paternal_in }

    HIFIASM.out.hap2_contigs
        .map { meta, fasta -> [ meta + [ 'haplotype': 2 ], fasta ] }
        .set { ch_gfastats_maternal_in }

    GFASTATS(
        ch_gfastats_paternal_in.mix(ch_gfastats_maternal_in),
        'fasta',
        '',
        '',
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(GFASTATS.out.versions)

    emit:
    assembled_haplotypes = GFASTATS.out.assembly // channel: [ val(meta), path(fasta) ]
    versions = ch_versions                       // channel: [ versions.yml ]
}
