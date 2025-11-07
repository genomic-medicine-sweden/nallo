
include { BCFTOOLS_MERGE        } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_SORT         } from '../../../modules/nf-core/bcftools/sort/main'
include { LONGPHASE_HAPLOTAG    } from '../../../modules/nf-core/longphase/haplotag/main'
include { LONGPHASE_PHASE       } from '../../../modules/nf-core/longphase/phase/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { SPLIT_MULTISAMPLE_VCF } from '../../../subworkflows/local/split_multisample_vcf/main'

workflow LONGPHASE {
    take:
    ch_snv_vcf       // channel: [ val(meta), path(vcf) ]
    ch_sv_vcf        // channel: [ val(meta), path(vcf) ] Optional
    ch_bam_bai       // channel: [ val(meta), path(bam), path(bai) ]
    fasta            // channel: [ val(meta), path(fasta) ]
    fai              // channel: [ val(meta), path(fai) ]
    phase_with_svs   // bool: Whether to include SVs in phasing (true) or not (false)

    main:
    ch_versions = Channel.empty()

    ch_snv_vcf
        .map { meta, vcf -> [ meta + [variant_type: 'snv'], vcf ] }
        .set { ch_snv_with_type }

    if (phase_with_svs) {
        ch_sv_vcf
            .map { meta, vcf -> [ meta + [variant_type: 'sv'], vcf ] }
            .mix( ch_snv_with_type )
            .set { ch_split_in }
    } else {
        ch_split_in = ch_snv_with_type
    }

    SPLIT_MULTISAMPLE_VCF (
        ch_split_in
    )
    ch_versions = ch_versions.mix(SPLIT_MULTISAMPLE_VCF.out.versions)

    SPLIT_MULTISAMPLE_VCF.out.split_vcf
        .branch { meta, vcf ->
            sv: meta.variant_type == 'sv'
                [meta - meta.subMap('variant_type'), vcf ]
            snv: meta.variant_type == 'snv'
                [meta - meta.subMap('variant_type'), vcf ]
        }
        .set { ch_split_vcfs }

    ch_bam_bai
        .map { meta, bam, bai -> [ [ id : meta.id, family_id : meta.family_id ], meta, bam, bai ] }
        .join( ch_split_vcfs.snv, failOnMismatch:true, failOnDuplicate:true )
        .set { ch_bam_vcf }

    if (phase_with_svs) {
        ch_bam_vcf
            .join( ch_split_vcfs.sv, failOnMismatch:true, failOnDuplicate:true )
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

    LONGPHASE_PHASE.out.snv_vcf
        .map { meta, vcf -> [ meta + [ variant_type: 'snv'], vcf ] }
        .mix ( LONGPHASE_PHASE.out.sv_vcf.map { meta, vcf -> [ meta + [ variant_type: 'sv'], vcf ] })
        .set { ch_bcftools_sort_in }

    // Sort all phased VCFs, ignoring variant types.
    BCFTOOLS_SORT( ch_bcftools_sort_in )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    BCFTOOLS_SORT.out.vcf
        .join (BCFTOOLS_SORT.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, vcf, tbi -> [ meta + [ id: meta.family_id ] - meta.subMap('family_id'), vcf, tbi ] }
        .groupTuple()
        .set { ch_phased_vcf }

    BCFTOOLS_MERGE (ch_phased_vcf, fasta, fai, [ [], [] ])
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    BCFTOOLS_MERGE.out.vcf
        .branch { meta, vcf ->
            snv: meta.variant_type == 'snv'
                [ meta - meta.subMap('variant_type'), vcf ]
            sv:  meta.variant_type == 'sv'
                [ meta - meta.subMap('variant_type'), vcf ]
        }
        .set { ch_phased_family_vcfs }


    BCFTOOLS_MERGE.out.index
        .branch { meta, tbi ->
            snv: meta.variant_type == 'snv'
                [ meta - meta.subMap('variant_type'), tbi ]
            sv:  meta.variant_type == 'sv'
                [ meta - meta.subMap('variant_type'), tbi ]
        }
        .set { ch_phased_family_vcf_index }

    ch_phased_family_snvs = ch_phased_family_vcfs.snv
    ch_phased_family_snvs_tbi = ch_phased_family_vcf_index.snv
    ch_phased_family_svs = ch_phased_family_vcfs.sv
    ch_phased_family_svs_tbi = ch_phased_family_vcf_index.sv


    // Haplotagging
    if (phase_with_svs) {
        LONGPHASE_PHASE.out.snv_vcf
            .join( LONGPHASE_PHASE.out.sv_vcf, failOnMismatch:true, failOnDuplicate:true )
            .map { meta, snvs, svs -> [ meta, snvs, svs, [] ] }
            .set { ch_vcfs_for_haplotag }
    } else {
        LONGPHASE_PHASE.out.snv_vcf
            .map { meta, vcf -> [ meta, vcf, [], [] ] }
            .set { ch_vcfs_for_haplotag }
    }

    LONGPHASE_HAPLOTAG (
        ch_bam_bai.join(ch_vcfs_for_haplotag, failOnMismatch:true, failOnDuplicate:true),
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(LONGPHASE_HAPLOTAG.out.versions)

    SAMTOOLS_INDEX (
        LONGPHASE_HAPLOTAG.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    LONGPHASE_HAPLOTAG.out.bam
        .join( SAMTOOLS_INDEX.out.bai, failOnMismatch:true, failOnDuplicate:true )
        .set { ch_bam_bai_haplotagged }

    emit:
    phased_family_snvs     = ch_phased_family_snvs      // channel: [ val(meta), path(vcf) ]
    phased_family_snvs_tbi = ch_phased_family_snvs_tbi  // channel: [ val(meta), path(tbi) ]
    phased_family_svs      = ch_phased_family_svs       // channel: [ val(meta), path(vcf) ]
    phased_family_svs_tbi  = ch_phased_family_svs_tbi   // channel: [ val(meta), path(tbi) ]
    haplotagged_bam_bai    = ch_bam_bai_haplotagged     // channel: [ val(meta), path(bam), path(bai) ]
    versions               = ch_versions                // channel: [ path(versions.yml) ]
}
