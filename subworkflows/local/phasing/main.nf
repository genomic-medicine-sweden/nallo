include { BCFTOOLS_CONCAT                                    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_LONGPHASE_SNV     } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_LONGPHASE_SV      } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_SORT                                      } from '../../../modules/nf-core/bcftools/sort/main'
include { CRAMINO as CRAMINO_PHASED                          } from '../../../modules/nf-core/cramino/main'
include { HIPHASE                                            } from '../../../modules/local/hiphase/main'
include { LONGPHASE_HAPLOTAG                                 } from '../../../modules/nf-core/longphase/haplotag/main'
include { LONGPHASE_PHASE                                    } from '../../../modules/nf-core/longphase/phase/main'
include { SAMTOOLS_CONVERT                                   } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_LONGPHASE         } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_WHATSHAP          } from '../../../modules/nf-core/samtools/index/main'
include { SPLIT_MULTISAMPLE_VCF as SPLIT_MULTISAMPLE_VCF_SNV } from '../../../subworkflows/local/split_multisample_vcf/main'
include { SPLIT_MULTISAMPLE_VCF as SPLIT_MULTISAMPLE_VCF_SV  } from '../../../subworkflows/local/split_multisample_vcf/main'
include { WHATSHAP_HAPLOTAG                                  } from '../../../modules/local/whatshap/haplotag/main'
include { WHATSHAP_PHASE                                     } from '../../../modules/local/whatshap/phase/main'
include { WHATSHAP_STATS                                     } from '../../../modules/local/whatshap/stats/main'

workflow PHASING {
    take:
    ch_snv_vcf       // channel: [ val(meta), path(vcf) ]
    ch_snv_vcf_index // channel: [ val(meta), path(tbi) ]
    ch_sv_vcf        // channel: [ val(meta), path(vcf) ] Optional
    ch_sv_vcf_index  // channel: [ val(meta), path(tbi) ] Optional
    ch_bam_bai       // channel: [ val(meta), path(bam), path(bai) ]
    fasta            // channel: [ val(meta), path(fasta) ]
    fai              // channel: [ val(meta), path(fai) ]
    phaser           //  string: Phasing tool to use
    phase_with_svs   //    bool: Whether to include SVs in phasing (true) or not (false)
    cram_output      //    bool: Publish alignments as CRAM (true) or BAM (false)

    main:
    ch_versions            = Channel.empty()

    // Phase variants and haplotag reads with Longphase
    if (phaser.equals("longphase")) {

        // The VCF channels are missing meta information for downstream subworkflows. So we preserve it during the join.
        SPLIT_MULTISAMPLE_VCF_SNV (
            ch_snv_vcf,
            "_snv"
        )
        ch_versions = ch_versions.mix(SPLIT_MULTISAMPLE_VCF_SNV.out.versions)

        ch_bam_bai
            .map { meta, bam, bai -> [ [ id : meta.id, family_id : meta.family_id ], meta, bam, bai ] }
            .join( SPLIT_MULTISAMPLE_VCF_SNV.out.split_vcf, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_bam_vcf }

        if (phase_with_svs) {

            SPLIT_MULTISAMPLE_VCF_SV (
                ch_sv_vcf,
                "_sv"
            )
            ch_versions = ch_versions.mix(SPLIT_MULTISAMPLE_VCF_SV.out.versions)

            ch_bam_vcf
                .join( SPLIT_MULTISAMPLE_VCF_SV.out.split_vcf, failOnMismatch:true, failOnDuplicate:true )
                .map { _meta, meta2, bam, bai, snvs, svs -> [ meta2, bam, bai, snvs, svs, [] ] }
                .set { ch_longphase_phase_in }

        } else {
            ch_bam_vcf
                .map { _meta, meta2, bam, bai, snvs -> [ meta2, bam, bai, snvs, [], [] ] }
                .set { ch_longphase_phase_in }
        }

        LONGPHASE_PHASE (
            ch_longphase_phase_in,
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)

        LONGPHASE_PHASE.out.vcf
            .transpose()
            .map { meta, vcf -> [ meta + [ sv : vcf.simpleName.endsWith("_SV") ], vcf ] }
            .set { ch_bcftools_sort_in }
        // Sort all phased VCFs, ignoring variant types.
        BCFTOOLS_SORT( ch_bcftools_sort_in )
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)


        // Separate phased SV and SNV VCFs and group by family
        BCFTOOLS_SORT.out.vcf
            .join (BCFTOOLS_SORT.out.tbi, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, vcf, tbi ->
                [ meta + [ id: meta.family_id ], vcf, tbi ]
            }
            .groupTuple()
            .branch { meta, _vcf, _tbi ->
                sv: meta.sv
                snv: !meta.sv
            }
            .set { ch_phased_vcf }

        BCFTOOLS_MERGE_LONGPHASE_SNV (ch_phased_vcf.snv, fasta, fai, [ [], [] ])
        ch_versions.mix(BCFTOOLS_MERGE_LONGPHASE_SNV.out.versions)

        ch_phased_family_snvs = BCFTOOLS_MERGE_LONGPHASE_SNV.out.vcf
        ch_phased_family_snvs_tbi = BCFTOOLS_MERGE_LONGPHASE_SNV.out.index

        // Although the following SV operations would be safe to run unconditionally,
        // we check whether we have any SVs to phase to avoid unnecessary concatenation
        if (phase_with_svs) {
            BCFTOOLS_MERGE_LONGPHASE_SV (ch_phased_vcf.sv, fasta, fai, [ [], [] ])
            ch_versions.mix(BCFTOOLS_MERGE_LONGPHASE_SV.out.versions)

            ch_phased_family_svs = BCFTOOLS_MERGE_LONGPHASE_SV.out.vcf
            ch_phased_family_svs_tbi = BCFTOOLS_MERGE_LONGPHASE_SV.out.index

            // Concatenate SNV and SV phased VCFs for Whatshap stats
            // If there are no phased SVs, we mix in an empty channel implicitly
            BCFTOOLS_MERGE_LONGPHASE_SNV.out.vcf
                .join(BCFTOOLS_MERGE_LONGPHASE_SNV.out.index, failOnMismatch: true, failOnDuplicate: true)
                .mix(BCFTOOLS_MERGE_LONGPHASE_SV.out.vcf
                    .join(BCFTOOLS_MERGE_LONGPHASE_SV.out.index, failOnMismatch: true, failOnDuplicate: true)
                )
                .map { meta, vcf, tbi -> [ meta - meta.subMap("sv"), vcf, tbi ] }
                .groupTuple()
                .set { ch_bcftools_concat_in }

            BCFTOOLS_CONCAT( ch_bcftools_concat_in )
            ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

            BCFTOOLS_CONCAT.out.vcf
                .join( BCFTOOLS_CONCAT.out.tbi )
                .set { ch_phased_vcf_index }
        } else {
            ch_phased_family_svs = ch_sv_vcf
            ch_phased_family_svs_tbi = ch_sv_vcf_index

            BCFTOOLS_MERGE_LONGPHASE_SNV.out.vcf
                .join(BCFTOOLS_MERGE_LONGPHASE_SNV.out.index, failOnMismatch: true, failOnDuplicate: true)
                .set { ch_phased_vcf_index }
        }

        // New if block here because this is concerned with haplotagging
        if (phase_with_svs) {
            LONGPHASE_PHASE.out.vcf
                .map { meta, vcfs ->
                    vcfs[1].simpleName.endsWith("_SV")
                        ? [ meta, vcfs[0], vcfs[1], [] ]
                        : [ meta, vcfs[1], vcfs[0], [] ]
                }
                .set { ch_vcfs_for_haplotag }
        } else {
            LONGPHASE_PHASE.out.vcf
                .map { meta, vcf -> [ meta, vcf, [], [] ] }
                .set { ch_vcfs_for_haplotag }
        }

        LONGPHASE_HAPLOTAG (
            ch_bam_bai.join(ch_vcfs_for_haplotag, failOnMismatch:true, failOnDuplicate:true),
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(LONGPHASE_HAPLOTAG.out.versions)

        SAMTOOLS_INDEX_LONGPHASE (
            LONGPHASE_HAPLOTAG.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_LONGPHASE.out.versions)

        LONGPHASE_HAPLOTAG.out.bam
            .join( SAMTOOLS_INDEX_LONGPHASE.out.bai, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_bam_bai_haplotagged }

    // Phase variants and haplotag reads with whatshap
    } else if (phaser.equals("whatshap")) {

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
            .set { ch_whatshap_out_split }

        ch_phased_family_snvs = ch_whatshap_out_split.vcf
        ch_phased_family_snvs_tbi = ch_whatshap_out_split.tbi
        ch_phased_family_svs = ch_sv_vcf
        ch_phased_family_svs_tbi = ch_sv_vcf_index

        ch_bam_bai
            .map { meta, bam, bai -> [ [ id : meta.family_id ], meta, bam, bai ] }
            .combine( WHATSHAP_PHASE.out.vcf_tbi, by: 0)
            .map { _meta, full_meta, bam, bai, vcf, tbi -> [ full_meta, vcf, tbi, bam, bai ] }
            .set { ch_whatshap_haplotag_in }

        WHATSHAP_HAPLOTAG (
            ch_whatshap_haplotag_in,
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)

        SAMTOOLS_INDEX_WHATSHAP (
            WHATSHAP_HAPLOTAG.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_WHATSHAP.out.versions)

        WHATSHAP_HAPLOTAG.out.bam
            .join( SAMTOOLS_INDEX_WHATSHAP.out.bai, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_bam_bai_haplotagged }

        WHATSHAP_PHASE.out.vcf_tbi
            .set { ch_phased_vcf_index }

    // Phase variants and haplotag reads with HiPhase
    } else if (phaser.equals("hiphase")) {

        ch_snv_vcf
            .join( ch_snv_vcf_index, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_snv_vcf_tbi }

        ch_bam_bai
            .map { meta, bam, bai -> [ [id : meta.family_id ], meta.id, bam, bai ]}
            .groupTuple()
            .map { meta, ids, bams, bais -> [ meta + [sample_ids: ids.toSet() ], bams, bais ]}
            .join( ch_snv_vcf_tbi, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_hiphase_bam_snv }

        if (phase_with_svs) {
            ch_sv_vcf
                .join( ch_sv_vcf_index, failOnMismatch:true, failOnDuplicate:true )
                .set { ch_sv_vcf_tbi }

            ch_hiphase_bam_snv
                .join( ch_sv_vcf_tbi, failOnMismatch: true, failOnDuplicate:true )
                .set { ch_hiphase_snv_in }
        } else {
            ch_hiphase_bam_snv
                .map { meta, bams, bais, snv_vcf, snv_tbi -> [ meta, bams, bais, snv_vcf, snv_tbi, [] , [] ] }
                .set { ch_hiphase_snv_in }
        }

        HIPHASE (
            ch_hiphase_snv_in,
            fasta,
            fai,
            true
        )
        ch_versions = ch_versions.mix(HIPHASE.out.versions)

        ch_phased_family_snvs = HIPHASE.out.vcfs
        ch_phased_family_snvs_tbi = HIPHASE.out.vcfs_tbi
        ch_phased_family_svs = phase_with_svs ? HIPHASE.out.sv_vcfs : ch_sv_vcf
        ch_phased_family_svs_tbi = phase_with_svs ? HIPHASE.out.sv_vcfs_tbi : ch_sv_vcf_index

        HIPHASE.out.bams
            .join( HIPHASE.out.bais, failOnMismatch:true, failOnDuplicate:true )
            .transpose()
            .combine(ch_bam_bai)
            .filter { _meta_phased, bam_phased, _bai_phased, meta_orig, _bam_orig, _bai_orig ->
                bam_phased.simpleName.startsWith(meta_orig.id)
            } // join does not allow arbitrary predicates, so we get the cross product and filter
            .map { _meta_phased, bam_phased, bai_phased, meta_orig, _bam_orig, _bai_orig ->
                [ meta_orig, bam_phased, bai_phased ]
            }
            .set { ch_bam_bai_haplotagged }

        HIPHASE.out.vcfs
            .join( HIPHASE.out.vcfs_tbi, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_phased_vcf_index }

        if (phase_with_svs) {
            HIPHASE.out.sv_vcfs
                .join( HIPHASE.out.sv_vcfs_tbi, failOnMismatch:true, failOnDuplicate:true )
                .mix( ch_phased_vcf_index )
                .groupTuple()
                .set { ch_bcftools_concat_in }

            BCFTOOLS_CONCAT (
                ch_bcftools_concat_in
            )
            ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
            BCFTOOLS_CONCAT.out.vcf
                .join( BCFTOOLS_CONCAT.out.tbi, failOnMismatch:true, failOnDuplicate:true )
                .set { ch_phased_vcf_index }
        }
    }

    // Phasing stats
    WHATSHAP_STATS ( ch_phased_vcf_index )
    ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions)

    if (cram_output) {
        SAMTOOLS_CONVERT (
            ch_bam_bai_haplotagged,
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)
    }

    // Phasing QC
    CRAMINO_PHASED ( ch_bam_bai_haplotagged )
    ch_versions = ch_versions.mix(CRAMINO_PHASED.out.versions)

    emit:
    phased_family_snvs     = ch_phased_family_snvs      // channel: [ val(meta), path(vcf) ]
    phased_family_snvs_tbi = ch_phased_family_snvs_tbi  // Channel: [ val(meta), path(tbi) ]
    phased_family_svs      = ch_phased_family_svs       // channel: [ val(meta), path(vcf) ]
    phased_family_svs_tbi  = ch_phased_family_svs_tbi   // Channel: [ val(meta), path(tbi) ]
    haplotagged_bam_bai    = ch_bam_bai_haplotagged     // channel: [ val(meta), path(bam), path(bai) ]
    stats                  = WHATSHAP_STATS.out.stats   // channel: [ val(meta), path(txt) ]
    versions               = ch_versions                // channel: [ path(versions.yml) ]
}
