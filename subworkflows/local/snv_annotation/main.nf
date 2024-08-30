include { ANNOTATE_CADD                               } from '../annotate_cadd/main'
include { ECHTVAR_ANNO                                } from '../../../modules/local/echtvar/anno/main'
include { BCFTOOLS_FILLTAGS as BCFTOOLS_FILLTAGS_ANNO } from '../../../modules/local/bcftools/filltags/main'
include { ENSEMBLVEP_VEP                              } from '../../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX as TABIX_VEP                    } from '../../../modules/nf-core/tabix/tabix/main'

workflow SNV_ANNOTATION {

    take:
    ch_vcf                // channel [mandatory] [ val(meta), path(vcf) ]
    ch_databases          // channel: [mandatory] [ val(meta), path(db) ]
    ch_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai                // channel: [mandatory] [ val(meta), path(fai) ]
    ch_vep_cache          // channel: [mandatory] [ path(cache) ]
    val_vep_cache_version // string: [mandatory] default: 110
    ch_vep_extra_files    // channel: [mandatory] [ path(files) ]
    val_annotate_cadd     // bool: [mandatory]
    ch_cadd_header        // channel: [mandatory] [ path(txt) ]
    ch_cadd_resources     // channel: [mandatory] [ path(annotation) ]
    ch_cadd_prescored     // channel: [mandatory] [ path(prescored) ]

    main:
    ch_versions = Channel.empty()
    ch_vep_in   = Channel.empty()

    // Annotate with chosen databases (GNOMAD,CADD + SAMPLES_DB)
    ECHTVAR_ANNO ( ch_vcf, ch_databases )
    ch_versions = ch_versions.mix(ECHTVAR_ANNO.out.versions)

    BCFTOOLS_FILLTAGS_ANNO(ECHTVAR_ANNO.out.bcf)
    ch_versions = ch_versions.mix(BCFTOOLS_FILLTAGS_ANNO.out.versions)

    // Annotating with CADD
    if (val_annotate_cadd) {
        ANNOTATE_CADD (
            ch_fai,
            BCFTOOLS_FILLTAGS_ANNO.out.vcf,
            BCFTOOLS_FILLTAGS_ANNO.out.tbi,
            ch_cadd_header,
            ch_cadd_resources,
            ch_cadd_prescored
        )
        ch_versions = ch_versions.mix(ANNOTATE_CADD.out.versions)

        ANNOTATE_CADD.out.vcf
            .map { meta, vcf -> [ meta, vcf, [] ] }
            .set { ch_vep_in }

    } else {
        BCFTOOLS_FILLTAGS_ANNO.out.vcf
            .map { meta, vcf -> [ meta, vcf, [] ] }
            .set { ch_vep_in }

    }

    ENSEMBLVEP_VEP (
        ch_vep_in,
        "GRCh38",
        "homo_sapiens",
        val_vep_cache_version,
        ch_vep_cache,
        ch_fasta,
        ch_vep_extra_files
    )
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)

    TABIX_VEP ( ENSEMBLVEP_VEP.out.vcf )
    ch_versions = ch_versions.mix(TABIX_VEP.out.versions)

    emit:
    vcf      = ENSEMBLVEP_VEP.out.vcf
    tbi      = TABIX_VEP.out.tbi
    versions = ch_versions
}
