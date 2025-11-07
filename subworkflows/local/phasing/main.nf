include { BCFTOOLS_CONCAT           } from '../../../modules/nf-core/bcftools/concat/main'
include { CRAMINO as CRAMINO_PHASED } from '../../../modules/nf-core/cramino/main'
include { HIPHASE                   } from '../../../modules/local/hiphase/main'
include { LONGPHASE                 } from '../../../subworkflows/local/longphase/main'
include { SAMTOOLS_CONVERT          } from '../../../modules/nf-core/samtools/convert/main'
include { WHATSHAP                  } from '../../../subworkflows/local/whatshap/main'
include { WHATSHAP_STATS            } from '../../../modules/local/whatshap/stats/main'

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

        LONGPHASE (
            ch_snv_vcf,
            ch_sv_vcf,
            ch_bam_bai,
            fasta,
            fai,
            phase_with_svs
        )
        ch_versions = ch_versions.mix(LONGPHASE.out.versions)

        ch_phased_family_snvs     = LONGPHASE.out.phased_family_snvs
        ch_phased_family_snvs_tbi = LONGPHASE.out.phased_family_snvs_tbi
        ch_phased_family_svs      = phase_with_svs ? LONGPHASE.out.phased_family_svs : ch_sv_vcf
        ch_phased_family_svs_tbi  = phase_with_svs ? LONGPHASE.out.phased_family_svs_tbi : ch_sv_vcf_index
        ch_bam_bai_haplotagged    = LONGPHASE.out.haplotagged_bam_bai

        // Set up concatenated VCF for stats
        if (phase_with_svs) {
            ch_phased_family_snvs
                .join(ch_phased_family_snvs_tbi, failOnMismatch: true, failOnDuplicate: true)
                .mix(ch_phased_family_svs
                    .join(ch_phased_family_svs_tbi, failOnMismatch: true, failOnDuplicate: true)
                )
                .groupTuple()
                .set { ch_bcftools_concat_in }

            BCFTOOLS_CONCAT( ch_bcftools_concat_in )
            ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

            BCFTOOLS_CONCAT.out.vcf
                .join( BCFTOOLS_CONCAT.out.tbi )
                .set { ch_phased_vcf_index }
        } else {
            ch_phased_family_snvs
                .join(ch_phased_family_snvs_tbi, failOnMismatch: true, failOnDuplicate: true)
                .set { ch_phased_vcf_index }
        }

    } else if (phaser.equals("whatshap")) {

        WHATSHAP(
            ch_snv_vcf,
            ch_bam_bai,
            fasta,
            fai
        )

        ch_versions = ch_versions.mix(WHATSHAP.out.versions)
        ch_phased_family_snvs     = WHATSHAP.out.phased_family_snvs
        ch_phased_family_snvs_tbi = WHATSHAP.out.phased_family_snvs_tbi
        ch_phased_family_svs      = ch_sv_vcf
        ch_phased_family_svs_tbi  = ch_sv_vcf_index
        ch_bam_bai_haplotagged    = WHATSHAP.out.haplotagged_bam_bai
        WHATSHAP.out.phased_family_snvs
            .join( WHATSHAP.out.phased_family_snvs_tbi, failOnMismatch:true, failOnDuplicate:true )
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
