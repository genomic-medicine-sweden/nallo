// TODO: BCFTOOLS processes should have unique names so that they are not used multiple times in other workflows?
include { ECHTVAR_ANNO                                  } from '../../modules/local/echtvar/anno/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SINGLESAMPLE   } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_SINGLESAMPLE } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_FILLTAGS as BCFTOOLS_FILLTAGS_ANNO   } from '../../modules/local/bcftools/filltags/main'
include { ENSEMBLVEP_VEP                                } from '../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX as TABIX_VEP                      } from '../../modules/nf-core/tabix/tabix/main'

workflow SNV_ANNOTATION {

    take:
    ch_bcf
    ch_single_sample_vcf
    ch_databases
    ch_fasta
    ch_vep_cache
    val_vep_cache_version // string: [mandatory] default: 110

    main:
    ch_versions       = Channel.empty()

    // Index and normalize single sample vcfs
    BCFTOOLS_INDEX_SINGLESAMPLE(ch_single_sample_vcf)

    BCFTOOLS_NORM_SINGLESAMPLE(
        ch_single_sample_vcf.join(BCFTOOLS_INDEX_SINGLESAMPLE.out.csi),
        ch_fasta
    )

    // Annotate with chosen databases (GNOMAD,CADD + SAMPLES_DB)

    ECHTVAR_ANNO ( BCFTOOLS_NORM_SINGLESAMPLE.out.vcf, ch_databases )
    BCFTOOLS_FILLTAGS_ANNO(ECHTVAR_ANNO.out.bcf)

    vep_in = BCFTOOLS_FILLTAGS_ANNO.out.vcf.map{ meta, vcf -> return [meta, vcf, []]}

    // Annotate with VEP as well

    ENSEMBLVEP_VEP(
        vep_in,
        "GRCh38",
        "homo_sapiens",
        val_vep_cache_version,
        ch_vep_cache,
        ch_fasta,
        []
    )

    TABIX_VEP ( ENSEMBLVEP_VEP.out.vcf )

    // Get versions
    ch_versions     = ch_versions.mix(BCFTOOLS_INDEX_SINGLESAMPLE.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_NORM_SINGLESAMPLE.out.versions)
    ch_versions     = ch_versions.mix(ECHTVAR_ANNO.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_FILLTAGS_ANNO.out.versions)
    ch_versions     = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
    ch_versions     = ch_versions.mix(TABIX_VEP.out.versions)

    emit:
    versions       = ch_versions
}
