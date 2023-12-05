// TODO: BCFTOOLS processes should have unique names so that they are not used multiple times in other workflows?
include { ECHTVAR_ANNO                                  } from '../../modules/local/echtvar/anno/main'
include { ECHTVAR_ENCODE                                } from '../../modules/local/echtvar/encode/main'
include { BCFTOOLS_NORM                                 } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SINGLESAMPLE   } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_INDEX                                } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_SINGLESAMPLE } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_FILLTAGS                             } from '../../modules/local/bcftools/filltags/main'
include { BCFTOOLS_FILLTAGS as BCFTOOLS_FILLTAGS_ANNO   } from '../../modules/local/bcftools/filltags/main'
include { ENSEMBLVEP_VEP                                } from '../../modules/nf-core/ensemblvep/vep/main'

workflow SNV_ANNOTATION {

    take:
    ch_bcf
    ch_single_sample_vcf
    ch_databases
    ch_fasta
    ch_vep_cache

    main:
    ch_versions       = Channel.empty()

    // Add allele count tag to mutlisample vcf
    BCFTOOLS_FILLTAGS(ch_bcf)
    // Index and normalize multisample vcf
    BCFTOOLS_INDEX(BCFTOOLS_FILLTAGS.out.vcf)
    BCFTOOLS_NORM(BCFTOOLS_FILLTAGS.out.vcf.join(BCFTOOLS_INDEX.out.csi), ch_fasta)

    // Index and normalize single sample vcfs
    BCFTOOLS_INDEX_SINGLESAMPLE(ch_single_sample_vcf)
    ss = ch_single_sample_vcf.join(BCFTOOLS_INDEX_SINGLESAMPLE.out.csi)
    BCFTOOLS_NORM_SINGLESAMPLE(ss, ch_fasta)

    // Make a cohort database using mutisample vcf
    ECHTVAR_ENCODE(BCFTOOLS_NORM.out.vcf)

    // combine input databases with cohort database
    db = ch_databases.concat(ECHTVAR_ENCODE.out.db.map{it[1]}).collect()

    // Annotate with chosen databases (GNOMAD,CADD + SAMPLES_DB)

    ECHTVAR_ANNO(BCFTOOLS_NORM_SINGLESAMPLE.out.vcf, db)
    BCFTOOLS_FILLTAGS_ANNO(ECHTVAR_ANNO.out.bcf)

    vep_in = BCFTOOLS_FILLTAGS_ANNO.out.vcf.map{ meta, vcf -> return [meta, vcf, []]}

    // Annotate with VEP as well

    ENSEMBLVEP_VEP(
        vep_in,
        "GRCh38",
        "homo_sapiens",
        "110",
        ch_vep_cache,
        ch_fasta,
        []
    )

    // Get versions
    ch_versions     = ch_versions.mix(BCFTOOLS_FILLTAGS.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_INDEX.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_NORM.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_INDEX_SINGLESAMPLE.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_NORM_SINGLESAMPLE.out.versions)
    ch_versions     = ch_versions.mix(ECHTVAR_ENCODE.out.versions)
    ch_versions     = ch_versions.mix(ECHTVAR_ANNO.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_FILLTAGS_ANNO.out.versions)
    ch_versions     = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)

    emit:
    versions       = ch_versions
}
