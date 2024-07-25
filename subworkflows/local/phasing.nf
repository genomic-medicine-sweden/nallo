include { BCFTOOLS_FILLFROMFASTA                       } from '../../modules/local/bcftools/fillfromfasta/main'
include { BCFTOOLS_REHEADER                            } from '../../modules/nf-core/bcftools/reheader/main'
include { CRAMINO as CRAMINO_PHASED                    } from '../../modules/local/cramino'
include { HIPHASE as HIPHASE_SNV                       } from '../../modules/local/hiphase/main'
include { HIPHASE as HIPHASE_SV                        } from '../../modules/local/hiphase/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_WHATSHAP    } from '../../modules/nf-core/samtools/index/main'
include { TABIX_BGZIPTABIX                             } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX                                  } from '../../modules/nf-core/tabix/tabix/main'
include { WHATSHAP_HAPLOTAG                            } from '../../modules/local/whatshap/haplotag/main'
include { WHATSHAP_PHASE                               } from '../../modules/local/whatshap/phase/main'
include { WHATSHAP_STATS                               } from '../../modules/local/whatshap/stats/main'

workflow PHASING {
    take:
        ch_vcf     // channel: [ val(meta), vcf ]
        ch_sv_vcf  // channel: [ val(meta), vcf ]
        ch_bam_bai // channel: [ val(meta), bam, bai ]
        fasta      // channel: [ val(meta), fasta ]
        fai        // channel: [ val(meta), fai ]

    main:
        ch_versions = Channel.empty()
        ch_bam_bai_haplotagged = Channel.empty()

        TABIX_TABIX(ch_vcf)
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

        if (params.phaser.equals("whatshap")) {

            WHATSHAP_PHASE( ch_vcf.join(ch_bam_bai), fasta, fai )
            ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions)

            WHATSHAP_STATS( WHATSHAP_PHASE.out.vcf_tbi )
            ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions)

            WHATSHAP_PHASE.out.vcf_tbi
                .join(ch_bam_bai)
                .set { ch_whatshap_haplotag_in }

            WHATSHAP_HAPLOTAG(ch_whatshap_haplotag_in, fasta, fai)
            ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)

            SAMTOOLS_INDEX_WHATSHAP( WHATSHAP_HAPLOTAG.out.bam )
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_WHATSHAP.out.versions)

            WHATSHAP_HAPLOTAG
                .out.bam
                .join(SAMTOOLS_INDEX_WHATSHAP.out.bai)
                .set { ch_bam_bai_haplotagged }

        } else if (params.phaser.equals("hiphase_snv")) {
            ch_vcf
                .join(TABIX_TABIX.out.csi)
                .join(ch_bam_bai)
                .set { ch_hiphase_snv_in }

            HIPHASE_SNV( ch_hiphase_snv_in, fasta, fai, true )
            ch_versions = ch_versions.mix(HIPHASE_SNV.out.versions)

            HIPHASE_SNV.out.bams
                .join(HIPHASE_SNV.out.bais)
                .set { ch_bam_bai_haplotagged }

        } else if (params.phaser.equals("hiphase_sv")) {
            // Sniffles specific...
            BCFTOOLS_REHEADER(
                ch_sv_vcf
                    .map { meta, vcf -> [meta, vcf, [], []] },
                [[],[]]
            )
            ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

            // Might be that newer versions of HiPhase ignores certain SVs
            // if BCFTOOLS_FILLFROMFASTA is not run, instead of craching
            BCFTOOLS_FILLFROMFASTA(BCFTOOLS_REHEADER.out.vcf, fasta)
            ch_versions = ch_versions.mix(BCFTOOLS_FILLFROMFASTA.out.versions)

            TABIX_BGZIPTABIX(BCFTOOLS_FILLFROMFASTA.out.vcf)
            ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

            TABIX_BGZIPTABIX.out.gz_tbi
                .map { meta, gz, tbi -> [ meta, gz ] }
                .set { ch_sv_vcf }

            TABIX_BGZIPTABIX.out.gz_tbi
                .map { meta, gz, tbi -> [ meta, tbi ] }
                .set { ch_sv_tbi }

            ch_vcf
                .concat(ch_sv_vcf)
                .groupTuple()
                .set { ch_hiphase_vcf }

            TABIX_TABIX.out.csi
                .concat(ch_sv_tbi)
                .groupTuple()
                .set { ch_hiphase_tbi }

            ch_hiphase_vcf
                .join(ch_hiphase_tbi)
                .join(ch_bam_bai)
                .set { ch_hiphase_in }

            HIPHASE_SV( ch_hiphase_in, fasta, fai, true )
            ch_versions = ch_versions.mix(HIPHASE_SV.out.versions)

            HIPHASE_SV.out.bams
                .join(HIPHASE_SV.out.bais)
                .set { ch_bam_bai_haplotagged }
        }

        CRAMINO_PHASED( ch_bam_bai_haplotagged )
        ch_versions = ch_versions.mix(CRAMINO_PHASED.out.versions)

    emit:
    haplotagged_bam_bai = ch_bam_bai_haplotagged // channel: [ val(meta), bam, bai ]
    versions            = ch_versions            // channel: [ versions.yml ]
}
