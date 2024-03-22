include { WHATSHAP_PHASE                            } from '../../modules/local/whatshap/phase/main'
include { WHATSHAP_STATS                            } from '../../modules/local/whatshap/stats/main'
include { WHATSHAP_HAPLOTAG                         } from '../../modules/local/whatshap/haplotag/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_WHATSHAP } from '../../modules/nf-core/samtools/index/main'
include { CRAMINO as CRAMINO_PHASED                 } from '../../modules/local/cramino'

include { TABIX_BGZIPTABIX                          } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX                               } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_REHEADER                         } from '../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_FILLFROMFASTA                    } from '../../modules/local/bcftools/fillfromfasta/main'
include { HIPHASE_SV                                } from '../../modules/local/hiphase/sv/main'
include { HIPHASE_SNV                               } from '../../modules/local/hiphase/snv/main'

workflow PHASING {
    take:
        ch_vcf // SNPs
        ch_sv_vcf
        ch_bam_bai
        fasta
        fai
    main:
        ch_versions = Channel.empty()

        // Sniffles specific
        BCFTOOLS_REHEADER(
            ch_sv_vcf
                .map { meta, vcf -> [meta, vcf, [], []] },
            [[],[]]
        )

        BCFTOOLS_FILLFROMFASTA(BCFTOOLS_REHEADER.out.vcf, fasta)
        TABIX_BGZIPTABIX(BCFTOOLS_FILLFROMFASTA.out.vcf)
        // Sniffles out

        TABIX_TABIX(ch_vcf)

        ch_vcf
            .join(TABIX_TABIX.out.tbi)
            .join(TABIX_BGZIPTABIX.out.gz_tbi)
            .join(ch_bam_bai)
            .set{ ch_hiphase_sv_in }

        ch_vcf
            .join(TABIX_TABIX.out.tbi)
            .join(ch_bam_bai)
            .set{ ch_hiphase_snv_in}

        HIPHASE_SV( ch_hiphase_sv_in, fasta, fai )

        HIPHASE_SNV( ch_hiphase_snv_in, fasta, fai )


        // Phase VCF
        WHATSHAP_PHASE ( ch_vcf.join(ch_bam_bai), fasta, fai )
        // Get phased stats
        WHATSHAP_STATS ( WHATSHAP_PHASE.out.vcf_tbi )

        WHATSHAP_PHASE.out.vcf_tbi
            .join(ch_bam_bai)
            .set{ ch_whatshap_haplotag_in }

        // Haplotag reads
        WHATSHAP_HAPLOTAG(ch_whatshap_haplotag_in, fasta, fai)

        // Index reads
        SAMTOOLS_INDEX_WHATSHAP ( WHATSHAP_HAPLOTAG.out.bam )

        // Combine haplotagged bams with bai
        WHATSHAP_HAPLOTAG
            .out.bam
            .join(SAMTOOLS_INDEX_WHATSHAP.out.bai)
            .set{ch_bam_bai_haplotagged}

         // Prepare inputs
        ch_bam_bai_haplotagged
            .map{ meta, bam, bai -> [ meta, bam, bai, [] ] }
            .set{ ch_mosdepth_in }

        CRAMINO_PHASED(ch_bam_bai_haplotagged)

        // Get versions
        ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions.first())
        ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions.first())
        ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_WHATSHAP.out.versions.first())
        ch_versions = ch_versions.mix(CRAMINO_PHASED.out.versions.first())

    emit:
    haplotagged_bam_bai = ch_bam_bai_haplotagged // channel: [ val(meta), bam, bai ]
    versions            = ch_versions            // channel: [ versions.yml ]
}
