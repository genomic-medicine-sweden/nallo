include { TRGT                                  } from '../../modules/local/trgt.nf'
include { TRVZ                                  } from '../../modules/local/trvz.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRGT } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRGT   } from '../../modules/nf-core/samtools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_TRGT   } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_TRGT } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_TRGT } from '../../modules/nf-core/bcftools/query/main'
include { GAWK                                  } from '../../modules/nf-core/gawk/main'

workflow REPEAT_ANALYSIS {

    take:
    ch_bam_bai
    ch_fasta
    ch_ped
    ch_trgt_bed

    main:
    ch_repeat_calls_vcf = Channel.empty()
    ch_versions     = Channel.empty()
    
    // This should work both with and without --trio, but really should skip PED-file and just rely on samplesheet
    ch_bam_bai
        .map{ [it[0]['id'], it[0], it[1], it[2] ]}
        .join(ch_ped.map{ [ it['id'], it['sex'] ] }, by: 0)
        .map{ [it[1], it[2], it[3], it[4]]}
        .set{ ch_tgrt_input } 
    
    TRGT ( ch_tgrt_input, ch_fasta, ch_trgt_bed )
    SAMTOOLS_SORT_TRGT(TRGT.out.bam)
    SAMTOOLS_INDEX_TRGT(SAMTOOLS_SORT_TRGT.out.bam)
   
    BCFTOOLS_SORT_TRGT(TRGT.out.vcf)
    BCFTOOLS_INDEX_TRGT(BCFTOOLS_SORT_TRGT.out.vcf)

    BCFTOOLS_SORT_TRGT.out.vcf
        .join(BCFTOOLS_INDEX_TRGT.out.csi)
        .set{ ch_bcftools_query_in }
    

    BCFTOOLS_QUERY_TRGT(ch_bcftools_query_in, [], [], [])
    
    GAWK(BCFTOOLS_QUERY_TRGT.out.txt, [])

    GAWK.out.output
        .transpose()
        .set{ ch_trvz_in }

    joiny = SAMTOOLS_SORT_TRGT.out.bam.join(SAMTOOLS_INDEX_TRGT.out.bai).join(ch_bcftools_query_in)
    //ch_trvz_in.cross(joiny).map{[it[0][0], it[1][1], it[1][2], it[1][3], it[1][4], it[0][1]]}.view()
    joiny
        .cross(ch_trvz_in)
        .map{[it[0][0], it[0][1], it[0][2], it[0][3], it[0][4], it[1][1]]}
        .set{ ch_trvz_in }

    TRVZ( ch_trvz_in, ch_fasta, ch_trgt_bed )

    

    ch_versions = ch_versions.mix(TRGT.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_TRGT.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TRGT.out.versions)
    
    emit:
    ch_repeat_calls_vcf        // channel: ? 
    versions = ch_versions // channel: [ versions.yml ]
}

