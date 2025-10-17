include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_LONGPHASE     } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_HIPHASE       } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_PLUGINSPLIT as BCFTOOLS_PLUGINSPLIT_SNV } from '../../../modules/nf-core/bcftools/pluginsplit/main'
include { BCFTOOLS_PLUGINSPLIT as BCFTOOLS_PLUGINSPLIT_SV  } from '../../../modules/nf-core/bcftools/pluginsplit/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_LONGPHASE       } from '../../../modules/nf-core/bcftools/merge/main'
include { CREATE_SPLIT_FILE as CREATE_SPLIT_FILE_SNV       } from '../../../modules/local/create_split_file/main'
include { CREATE_SPLIT_FILE as CREATE_SPLIT_FILE_SV        } from '../../../modules/local/create_split_file/main'
include { CRAMINO as CRAMINO_PHASED                        } from '../../../modules/local/cramino/main'
include { HIPHASE                                          } from '../../../modules/local/hiphase/main'
include { LONGPHASE_HAPLOTAG                               } from '../../../modules/nf-core/longphase/haplotag/main'
include { LONGPHASE_PHASE                                  } from '../../../modules/nf-core/longphase/phase/main'
include { SAMTOOLS_CONVERT                                 } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_LONGPHASE       } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_WHATSHAP        } from '../../../modules/nf-core/samtools/index/main'
include { WHATSHAP_HAPLOTAG                                } from '../../../modules/local/whatshap/haplotag/main'
include { WHATSHAP_PHASE                                   } from '../../../modules/local/whatshap/phase/main'
include { WHATSHAP_STATS                                   } from '../../../modules/local/whatshap/stats/main'

workflow PHASING {
    take:
    ch_snv_vcf       // channel: [ val(meta), path(vcf) ]
    ch_snv_vcf_index // channel: [ val(meta), path(tbi) ]
    ch_sv_vcf        // channel: [ val(meta), path(vcf) ]
    ch_sv_vcf_index  // channel: [ val(meta), path(tbi) ]
    ch_bam_bai       // channel: [ val(meta), path(bam), path(bai) ]
    fasta            // channel: [ val(meta), path(fasta) ]
    fai              // channel: [ val(meta), path(fai) ]
    phaser           //  string: Phasing tool to use
    cram_output      //    bool: Publish alignments as CRAM (true) or BAM (false)

    main:
    ch_versions            = Channel.empty()

    // Phase variants and haplotag reads with Longphase
    if (phaser.equals("longphase")) {

        ch_snv_vcf
            .map { meta, _vcf -> [ meta, meta.sample_ids ] }
            .set { ch_family_samples }

        // Create sample files for bcftools +split
        CREATE_SPLIT_FILE_SNV ( ch_family_samples, "_snv" )
        CREATE_SPLIT_FILE_SV ( ch_family_samples, "_sv" )

        // Collect the files' content into channels
        // We must merge the sample IDs based on file names later
        // And this allows us to forego any assumptions about the filenames
        CREATE_SPLIT_FILE_SNV.out.txt
            .splitCsv(sep: "\t", elem: 1)
            .set { ch_snv_split_names }

        CREATE_SPLIT_FILE_SV.out.txt
            .splitCsv(sep: "\t", elem: 1)
            .set { ch_sv_split_names }

        ch_snv_vcf
            .join(ch_snv_vcf_index, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_snv }

        ch_sv_vcf
            .join(ch_sv_vcf_index, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_sv }


        BCFTOOLS_PLUGINSPLIT_SNV (
            ch_snv,
            CREATE_SPLIT_FILE_SNV.out.txt,
            [],
            [],
            []
        )
        BCFTOOLS_PLUGINSPLIT_SNV.out.vcf
            .transpose()
            .map { meta, file -> [ file.simpleName, meta, file ]}
            .join(ch_snv_split_names, by: [0, 1], failOnMismatch: true, failOnDuplicate: true)
            .map { _filename, meta, file, sample -> [ meta - meta.subMap("id") + [ id: sample, family_id: meta.id ], file ] }
            .set { ch_split_snv_vcf }

        BCFTOOLS_PLUGINSPLIT_SV (
            ch_sv,
            CREATE_SPLIT_FILE_SV.out.txt,
            [],
            [],
            []
        )
        BCFTOOLS_PLUGINSPLIT_SV.out.vcf
            .transpose()
            .map { meta, file -> [ file.simpleName, meta, file ]}
            .join(ch_sv_split_names, by: [0, 1], failOnMismatch: true, failOnDuplicate: true)
            .map { _filename, meta, file, sample -> [ meta - meta.subMap("id") + [ id: sample, family_id: meta.id ], file ] }
            .set { ch_split_sv_vcf }

        ch_bam_bai
            .map { meta, bam, bai -> [ [ id : meta.id, family_id : meta.family_id ], meta, bam, bai ] }
            .join( ch_split_snv_vcf, failOnMismatch:true, failOnDuplicate:true )
            .join( ch_split_sv_vcf , remainder:true, failOnDuplicate:true ) // Will set svs to null in case there are none
            .map { _meta, meta2, bam, bai, snvs, svs -> [ meta2, bam, bai, snvs, svs ?: [], [] ] }
            .set { ch_longphase_phase_in }

        LONGPHASE_PHASE (
            ch_longphase_phase_in,
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)

        LONGPHASE_PHASE.out.vcf
            .map { meta, vcfs -> [ meta, vcfs, [] ]}
            .set { ch_bcftools_concat_in }

        BCFTOOLS_CONCAT_LONGPHASE( ch_bcftools_concat_in )
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_LONGPHASE.out.versions)

        BCFTOOLS_CONCAT_LONGPHASE.out.vcf
            .map { meta, vcfs -> [ [ id : meta.family_id ], meta.id, vcfs, ] }
            .groupTuple()
            .map { meta, ids, vcfs -> [ meta + [ sample_ids : ids ], vcfs.flatten() ]}
            .set { ch_bcftools_merge_in }

        // Longphase emits 2 VCFs if we supplied svs
        // Concat all VCFs for each sample for publishing and stats
        BCFTOOLS_MERGE_LONGPHASE ( ch_bcftools_merge_in, [], [], [] )
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE_LONGPHASE.out.versions)

        BCFTOOLS_MERGE_LONGPHASE.out.vcf
            .join(BCFTOOLS_MERGE_LONGPHASE.out.index, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_phased_vcf_index }

        // However, haplotag requires the VCFs to stay separate
        // We need to "flatten" the VCF list into separate values in the output tuple if we have SVs
        // We can identify which VCF is which by the file name. We should not rely on order in the list
        LONGPHASE_PHASE.out.vcf
            .map { meta, vcfs ->
                vcfs instanceof List
                    ? vcfs[1].baseName.endsWith("_SV")
                        ? [ meta, vcfs[0], vcfs[1] ]
                        : [ meta, vcfs[1], vcfs[0] ]
                    : [ meta, vcfs, [] ]
            }
            .set { ch_phased_vcf }

        ch_bam_bai
            .join( ch_phased_vcf, failOnMismatch:true, failOnDuplicate:true )
            .map { meta, bam, bai, snvs, svs -> [ meta, bam, bai, snvs, svs, [] ] }
            .set { ch_longphase_haplotag_in }

        LONGPHASE_HAPLOTAG (
            ch_longphase_haplotag_in,
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

        ch_bam_bai
            .map { meta, bam, bai -> [ [id : meta.family_id ], meta.id, bam, bai ]}
            .groupTuple()
            .map { meta, ids, bams, bais -> [ meta + [sample_ids: ids.toSet() ], bams, bais ]}
            .join( ch_snv_vcf, failOnMismatch:true, failOnDuplicate:true )
            .join( ch_snv_vcf_index, failOnMismatch:true, failOnDuplicate:true )
            .join( ch_sv_vcf, failOnMismatch:true, failOnDuplicate:true )
            .join( ch_sv_vcf_index, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_hiphase_snv_in }

        HIPHASE (
            ch_hiphase_snv_in,
            fasta,
            fai,
            true
        )
        ch_versions = ch_versions.mix(HIPHASE.out.versions)

        HIPHASE.out.bams
            .join( HIPHASE.out.bais, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_bam_bai_haplotagged }

        HIPHASE.out.vcfs
            .join( HIPHASE.out.vcfs_tbi, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_phased_snv_index }

        HIPHASE.out.sv_vcfs
            .join( HIPHASE.out.sv_vcfs_tbi, failOnMismatch:true, failOnDuplicate:true )
            .mix( ch_phased_snv_index )
            .groupTuple()
            .set { ch_bcftools_concat_in }

        BCFTOOLS_CONCAT_HIPHASE (
            ch_bcftools_concat_in
        )

        BCFTOOLS_CONCAT_HIPHASE.out.vcf
            .join( BCFTOOLS_CONCAT_HIPHASE.out.tbi, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_phased_vcf_index }

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
    haplotagged_bam_bai = ch_bam_bai_haplotagged   // channel: [ val(meta), path(bam), path(bai) ]
    stats               = WHATSHAP_STATS.out.stats // channel: [ val(meta), path(txt) ]
    versions            = ch_versions              // channel: [ path(versions.yml) ]
}
