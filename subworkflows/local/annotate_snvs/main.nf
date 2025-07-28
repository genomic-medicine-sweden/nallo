include { ANNOTATE_CADD                               } from '../annotate_cadd/main'
include { BCFTOOLS_FILLTAGS as BCFTOOLS_FILLTAGS_ANNO } from '../../../modules/local/bcftools/filltags/main'
include { BCFTOOLS_VIEW                               } from '../../../modules/nf-core/bcftools/view/main'
include { ECHTVAR_ANNO                                } from '../../../modules/local/echtvar/anno/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_SNV            } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX as TABIX_ENSEMBLVEP_SNV         } from '../../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_SNVS {

    take:
    ch_vcf                    // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_databases              // channel:  [optional] [ path(db) ]
    ch_fasta                  // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai                    // channel: [mandatory] [ val(meta), path(fai) ]
    ch_vep_cache              // channel: [mandatory] [ path(cache) ]
    val_vep_cache_version     //  string: [mandatory] default: 110
    ch_vep_extra_files        // channel: [mandatory] [ path(files) ]
    annotate_cadd             //    bool: [mandatory] should CADD be used to annotate indels
    annotate_echtvar          //    bool: [mandatory] should echtvar be used to annotate variants
    ch_cadd_header            // channel:  [optional] [ path(txt) ]
    ch_cadd_resources         // channel:  [optional] [ val(meta), path(annotation) ]
    ch_cadd_prescored_indels  // channel:  [optional] [ val(meta), path(prescored) ]
    pre_vep_filter            //    bool: [mandatory] should filtering be done before annotating with CADD and VEP

    main:
    ch_versions = Channel.empty()

    // Annotate with chosen databases
    if (annotate_echtvar) {
        ECHTVAR_ANNO (
            ch_vcf,
            ch_databases
        )
        ch_versions = ch_versions.mix(ECHTVAR_ANNO.out.versions)
    }

    // Allows for filtering before annotating with VEP
    if (pre_vep_filter) {
        (annotate_echtvar ? ECHTVAR_ANNO.out.bcf : ch_vcf)
            .map { meta, vcf -> [ meta, vcf, [] ] }
            .set { ch_bcftools_view_input }

        BCFTOOLS_VIEW (
            ch_bcftools_view_input,
            [],
            [],
            []
        )
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
    }

    //
    // Annotating with CADD
    //
    // If we don't have an index the subworkflow will create one, but needs to join on meta to create a [ meta, vcf, [] ] channel
    if (annotate_cadd) {
        ANNOTATE_CADD (
            ch_fai,
            pre_vep_filter ? BCFTOOLS_VIEW.out.vcf : annotate_echtvar ? ECHTVAR_ANNO.out.bcf : ch_vcf,
            pre_vep_filter ? BCFTOOLS_VIEW.out.tbi : ch_vcf.map { meta, _vcf -> [ meta, [] ] },
            ch_cadd_header,
            ch_cadd_resources,
            ch_cadd_prescored_indels
        )
        ch_versions = ch_versions.mix(ANNOTATE_CADD.out.versions)
    }

    (annotate_cadd ? ANNOTATE_CADD.out.vcf : pre_vep_filter ? BCFTOOLS_VIEW.out.vcf : annotate_echtvar ? ECHTVAR_ANNO.out.bcf : ch_vcf)
        .map { meta, vcf -> [ meta, vcf, [] ] }
        .set { ch_vep_in }

    // Always annotate with VEP
    ENSEMBLVEP_SNV (
        ch_vep_in,
        "GRCh38",
        "homo_sapiens",
        val_vep_cache_version,
        ch_vep_cache,
        ch_fasta,
        ch_vep_extra_files
    )
    ch_versions = ch_versions.mix(ENSEMBLVEP_SNV.out.versions)

    TABIX_ENSEMBLVEP_SNV (
        ENSEMBLVEP_SNV.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_ENSEMBLVEP_SNV.out.versions)

    emit:
    vcf      = ENSEMBLVEP_SNV.out.vcf
    tbi      = TABIX_ENSEMBLVEP_SNV.out.tbi
    versions = ch_versions
}
