include { MODKIT_SAMPLEPROBS                      } from '../../../modules/local/modkit/sample_probs/main'
include { SCATTER_GENOME as SCATTER_GENOME_METHYL } from '../../../subworkflows/local/scatter_genome'
include { MODKIT_PILEUP                           } from '../../../modules/local/modkit/pileup/main'
include { BEDMETHYL_CONCAT                        } from '../../../modules/local/bedmethyl_concat/main'
include { TABIX_BGZIPTABIX                        } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { MODKIT_BEDMETHYLTOBIGWIG                } from '../../../modules/nf-core/modkit/bedmethyltobigwig/main'

workflow METHYLATION {

    take:
    ch_bam_bai             // channel: [ val(meta), bam, bai ]
    ch_fasta               // channel: [ val(meta), fasta ]
    ch_fai                 // channel: [ val(meta), fai ]
    ch_bed                 // channel: [ val(meta), bed ]
    modcodes               // String or List

    main:
    ch_versions = Channel.empty()

    MODKIT_SAMPLEPROBS(ch_bam_bai)

    // Scatter genome for methyl too (same nb chunks/processes as for 'snv_calling')
    SCATTER_GENOME_METHYL (
        ch_fai,
        ch_bed,      // BED file to scatter
        !params.methylation_call_regions, // Make bed from fai
        !params.skip_methylation_pileups,
        params.snv_calling_processes
    )

    // Combine the BED intervals with BAM/BAI files to create a region-bam-bai for each sample.
    // This uses the whole BAM files for each region instead of splitting them.
    ch_bam_bai
        .join(MODKIT_SAMPLEPROBS.out.probs)
        .combine(SCATTER_GENOME_METHYL.out.bed_intervals)
        .map { meta, bam, bai, probs, bed, intervals ->
            [ meta + [ num_intervals: intervals, region: bed, probs: probs ], bam, bai, bed ]
        }
        .set { call_methyl_input }

    // Performs pileups per haplotype if the phasing workflow is on, set in config
    MODKIT_PILEUP (call_methyl_input, ch_fasta)
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)

    if(!params.skip_phasing) {
        MODKIT_PILEUP.out.modkit_H0
            .mix(MODKIT_PILEUP.out.modkit_H1, MODKIT_PILEUP.out.modkit_H2)
            .map{ meta, group, bedmethyl -> [meta - meta.subMap('region') + [ group: group ], bedmethyl] }
            .groupTuple(by: 0)
            .set { bedmethyl_to_concat_per_sample }
    } else {
        MODKIT_PILEUP.out.modkit_Hu
            .map{ meta, group, bedmethyl -> [meta - meta.subMap('region') + [ group: group ], bedmethyl] }
            .groupTuple(by: 0)
            .set { bedmethyl_to_concat_per_sample }
    }

    BEDMETHYL_CONCAT(bedmethyl_to_concat_per_sample)

    BEDMETHYL_CONCAT.out.sorted
        .transpose()
        .set { ch_bedmethyl }

    TABIX_BGZIPTABIX ( ch_bedmethyl )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    // Only convert files with content
    ch_bedmethyl
        .filter { _meta, bed -> bed.size() > 0 }
        .set { ch_bedmethyl_to_bigwig_in }

    MODKIT_BEDMETHYLTOBIGWIG ( ch_bedmethyl_to_bigwig_in, ch_fai, modcodes)
    ch_versions = ch_versions.mix(MODKIT_BEDMETHYLTOBIGWIG.out.versions)

    emit:
    bed      = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, bed, _tbi -> [ meta, bed ] } // channel: [ val(meta), path(bed) ]
    tbi      = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, _bed, tbi -> [ meta, tbi ] } // channel: [ val(meta), path(tbi) ]
    bigwig   = MODKIT_BEDMETHYLTOBIGWIG.out.bw
    versions = ch_versions // channel: [ versions.yml ]
}
