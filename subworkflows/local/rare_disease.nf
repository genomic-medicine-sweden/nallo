include { ECHTVAR_ANNO              } from '../../modules/local/echtvar/anno/main'
include { ECHTVAR_ENCODE              } from '../../modules/local/echtvar/encode/main'
include { BCFTOOLS_NORM             } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SINGLESAMPLE            } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_INDEX            } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_SINGLESAMPLE            } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_VIEW  } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_FILLTAGS  } from '../../modules/local/bcftools/filltags/main'
include { BCFTOOLS_FILLTAGS as BCFTOOLS_FILLTAGS_ANNO  } from '../../modules/local/bcftools/filltags/main'

workflow RARE_DISEASE {

    take:
    ch_bcf
    ch_single_sample_vcf
    ch_databases
    ch_fasta

    main:
    ch_versions       = Channel.empty()

    BCFTOOLS_FILLTAGS(ch_bcf)
    BCFTOOLS_INDEX(BCFTOOLS_FILLTAGS.out.vcf)
    BCFTOOLS_NORM(BCFTOOLS_FILLTAGS.out.vcf.join(BCFTOOLS_INDEX.out.csi), ch_fasta)

    BCFTOOLS_INDEX_SINGLESAMPLE(ch_single_sample_vcf)
    ss = ch_single_sample_vcf.join(BCFTOOLS_INDEX_SINGLESAMPLE.out.csi)
    BCFTOOLS_NORM_SINGLESAMPLE(ss, ch_fasta)

    // Make a cohort database
    ECHTVAR_ENCODE(BCFTOOLS_NORM.out.vcf)

    // combine input databases with cohort database
    db = ch_databases.concat(ECHTVAR_ENCODE.out.db.map{it[1]}).collect()

    //CADD

    ECHTVAR_ANNO(BCFTOOLS_NORM_SINGLESAMPLE.out.vcf, db)
    BCFTOOLS_FILLTAGS_ANNO(ECHTVAR_ANNO.out.bcf)


    // Get versions
    ch_versions     = ch_versions.mix(BCFTOOLS_FILLTAGS.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_INDEX.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_NORM.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_INDEX_SINGLESAMPLE.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_NORM_SINGLESAMPLE.out.versions)
    ch_versions     = ch_versions.mix(ECHTVAR_ENCODE.out.versions)
    ch_versions     = ch_versions.mix(ECHTVAR_ANNO.out.versions)
    ch_versions     = ch_versions.mix(BCFTOOLS_FILLTAGS_ANNO.out.versions)

    emit:
    versions       = ch_versions
}
