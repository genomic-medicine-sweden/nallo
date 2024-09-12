//
// Workflow to call and merge SNVs
//
include { BCFTOOLS_CONCAT                             } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_FILLTAGS                           } from '../../../modules/local/bcftools/filltags/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MULTISAMPLE  } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SINGLESAMPLE } from '../../../modules/nf-core/bcftools/norm/main'
include { DEEPVARIANT                                 } from '../../../modules/nf-core/deepvariant/main'
include { GLNEXUS                                     } from '../../../modules/nf-core/glnexus/main'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_bam_bai_bed // channel: [mandatory] [ val(meta), path(bam), path(bai), path(call_legion_bed) ]
    ch_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai         // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed         // channel: [optional] [ val(meta), path(input_bed) ]
    ch_par_bed     // channel: [mandatory] [ val(meta), path(par_bed) ]

    main:
    ch_versions = Channel.empty()

    ch_bam_bai_bed
        // Add call legion to meta so we can group by it later
        .map { meta, bam, bai, bed ->
            [ meta + [ 'legion': bed ], bam, bai, bed ]
        }
        .set { ch_deepvariant_in }

    DEEPVARIANT ( ch_deepvariant_in, ch_fasta, ch_fai, [[],[]], ch_par_bed )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

    // First remove legion so we can group per sample
    // Then after grouping remove num_intervals since to match the meta of other workflows
    DEEPVARIANT.out.vcf
        .map { meta, vcf ->
            new_meta = meta - meta.subMap('legion')
            [ groupKey(new_meta, new_meta.num_intervals ), vcf ]
        }
        .groupTuple()
        .join( DEEPVARIANT.out.vcf_tbi
            .map{ meta, tbi ->
                new_meta = meta - meta.subMap('legion')
                [ groupKey(new_meta, new_meta.num_intervals ), tbi ]
            }
            .groupTuple()
        )
        .map { meta, vcf, tbi ->
            [ meta - meta.subMap('num_intervals'), vcf, tbi ]
        }
        .set{ ch_concat_singlesample_in }

    // This creates a singlesample VCF containing ALL legions
    BCFTOOLS_CONCAT ( ch_concat_singlesample_in )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    // Which is then normalized, and ready to be used
    // in processes that require SNVs, but not annotated SNVs
    BCFTOOLS_NORM_SINGLESAMPLE ( BCFTOOLS_CONCAT.out.vcf.map { meta, vcf -> [ meta, vcf, [] ] }, ch_fasta )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM_SINGLESAMPLE.out.versions)

    // This creates a multisample VCF, with legions from ONE bed file
    DEEPVARIANT.out.gvcf
        .map { meta, gvcf ->
            [ meta.legion.name, meta.project, meta.phenotype == 2, gvcf ]
        }
        .groupTuple() // Group all files together per legion
        // If any of the samples in the VCF have an affected phenotype (2)
        // add this to the meta of the multisample VCF to know if we should run RANK_VARIANTS or not
        .map { meta, project, affected, gvcfs ->
            new_meta = [
                'id': meta,
                'project': project.first(), // Works only because only one project per run is allowed
                'contains_affected': affected.any(),
            ]
            [ new_meta, gvcfs ]
        }
        .set{ glnexus_in }

    GLNEXUS( glnexus_in, ch_bed )
    ch_versions = ch_versions.mix(GLNEXUS.out.versions)

    // Add allele count tag to multisample bcf
    BCFTOOLS_FILLTAGS ( GLNEXUS.out.bcf )
    ch_versions = ch_versions.mix(BCFTOOLS_FILLTAGS.out.versions)

    BCFTOOLS_FILLTAGS.out.vcf
        .map { meta, vcf -> [ meta, vcf, [] ] }
        .set { bcftools_norm_in }

    // Decompose and normalize variants
    BCFTOOLS_NORM_MULTISAMPLE ( bcftools_norm_in, ch_fasta )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM_MULTISAMPLE.out.versions)

    emit:
    snp_calls_vcf = BCFTOOLS_NORM_SINGLESAMPLE.out.vcf // channel: [ val(meta), path(bcf) ]
    combined_bcf  = BCFTOOLS_NORM_MULTISAMPLE.out.vcf  // channel: [ val(meta), path(bcf) ]
    combined_csi  = BCFTOOLS_NORM_MULTISAMPLE.out.csi  // channel: [ val(meta), path(csi) ]
    versions      = ch_versions                        // channel: [ path(versions.yml) ]
}
