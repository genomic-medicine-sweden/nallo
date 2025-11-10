include { HIPHASE } from '../../../modules/local/hiphase/main'

workflow RUN_HIPHASE {
    take:
    ch_snv_vcf       // channel: [ val(meta), path(vcf) ]
    ch_snv_vcf_index // channel: [ val(meta), path(tbi) ]
    ch_sv_vcf        // channel: [ val(meta), path(vcf) ] Optional
    ch_sv_vcf_index  // channel: [ val(meta), path(tbi) ] Optional
    ch_bam_bai       // channel: [ val(meta), path(bam), path(bai) ]
    fasta            // channel: [ val(meta), path(fasta) ]
    fai              // channel: [ val(meta), path(fai) ]
    phase_with_svs   // bool: Whether to include SVs in phasing (true) or not (false)

    main:
    ch_versions = Channel.empty()

    // Prepare SNV VCF with index
    ch_snv_vcf.dump(tag: 'SNV VCFs')
        .join( ch_snv_vcf_index, failOnMismatch:true, failOnDuplicate:true )
        .set { ch_snv_vcf_tbi }

    // Group BAM files by family and join with SNV VCF
    ch_bam_bai
        .map { meta, bam, bai -> [ [id : meta.family_id ], meta.id, bam, bai ]}
        .groupTuple()
        .map { meta, ids, bams, bais -> [ meta + [sample_ids: ids.toSet() ], bams, bais ]}.dump(tag:"Grouped BAMs")
        .join( ch_snv_vcf_tbi.dump(tag: 'VCFs'), failOnMismatch:true, failOnDuplicate:true )
        .set { ch_hiphase_bam_snv }

    // Prepare input based on whether SVs are included
    if (phase_with_svs) {
        ch_sv_vcf
            .join( ch_sv_vcf_index, failOnMismatch:true, failOnDuplicate:true )
            .set { ch_sv_vcf_tbi }

        ch_hiphase_bam_snv
            .join( ch_sv_vcf_tbi, failOnMismatch: true, failOnDuplicate:true )
            .set { ch_hiphase_input }
    } else {
        ch_hiphase_bam_snv
            .map { meta, bams, bais, snv_vcf, snv_tbi -> [ meta, bams, bais, snv_vcf, snv_tbi, [] , [] ] }
            .set { ch_hiphase_input }
    }

    // Run HiPhase
    HIPHASE (
        ch_hiphase_input,
        fasta,
        fai,
        true
    )
    ch_versions = ch_versions.mix(HIPHASE.out.versions)

    // Prepare haplotagged BAM output by matching with original metadata
    HIPHASE.out.bams
        .join( HIPHASE.out.bais, failOnMismatch:true, failOnDuplicate:true )
        .transpose()
        .combine(ch_bam_bai)
        .filter { _meta_phased, bam_phased, _bai_phased, meta_orig, _bam_orig, _bai_orig ->
            bam_phased.simpleName.startsWith(meta_orig.id)
        }
        .map { _meta_phased, bam_phased, bai_phased, meta_orig, _bam_orig, _bai_orig ->
            [ meta_orig, bam_phased, bai_phased ]
        }
        .set { ch_haplotagged_bam_bai }

    emit:
    phased_snvs             = HIPHASE.out.vcfs                                           // channel: [ val(meta), path(vcf) ]
    phased_snvs_tbi         = HIPHASE.out.vcfs_tbi                                       // channel: [ val(meta), path(tbi) ]
    phased_svs              = phase_with_svs ? HIPHASE.out.sv_vcfs : ch_sv_vcf           // channel: [ val(meta), path(vcf) ]
    phased_svs_tbi          = phase_with_svs ? HIPHASE.out.sv_vcfs_tbi : ch_sv_vcf_index // channel: [ val(meta), path(tbi) ]
    haplotagged_bam_bai     = ch_haplotagged_bam_bai                                          // channel: [ val(meta), path(bam), path(bai) ]
    versions                = ch_versions                                                     // channel: [ path(versions.yml) ]
}
