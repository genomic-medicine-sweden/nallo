include { DEEPVARIANT                      } from '../../modules/local/google/deepvariant'
include { PEPPER_MARGIN_DEEPVARIANT        } from '../../modules/local/pepper_margin_deepvariant'
include { DEEPTRIO                         } from '../../modules/local/google/deeptrio'
include { GLNEXUS                          } from '../../modules/nf-core/glnexus'
include { BCFTOOLS_VIEW_REGIONS            } from '../../modules/local/bcftools/view_regions'
include { TABIX_TABIX as TABIX_EXTRA_GVCFS } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DV } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_SINGLESAMPLE } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_CONCAT_SINGLESAMPLE } from '../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_DV } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_DV } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX as TABIX_CONCAT_UNMERGED_MULTISAMPLE } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_CONCAT_UNMERGED_SINGLESAMPLE } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_CONCAT_MERGED_SINGLESAMPLE } from '../../modules/nf-core/tabix/tabix/main'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_bam_bai
    ch_extra_gvcfs
    ch_fasta
    ch_fai
    ch_bed

    main:
    ch_snp_calls_vcf  = Channel.empty()
    ch_snp_calls_gvcf = Channel.empty()
    ch_combined_bcf   = Channel.empty()
    ch_versions       = Channel.empty()

    // Adding regions to ch_bam_bai surely messes up this below, therefore just remove it for now..?

    if (params.variant_caller == 'deeptrio') {
        ch_reads
            .map{it[0]}
            .combine(ch_bam_bai
                .map{ meta, bam, bai, region -> [meta.id, bam, bai]
                }
            )
            .branch{
                kid: it[0]['id']          == it[1]
                mom: it[0]['maternal_id'] == it[1]
                dad: it[0]['paternal_id'] == it[1]
            }
            .set{branch_result}

        branch_result
            .kid
            .join(branch_result.dad, remainder: true)
            .join(branch_result.mom, remainder: true)
            .map{[it[0], it[2], it[3], it[5], it[6], it[8], it[9]]} // [meta, kid_bam, kid_bai, dad_bam, dad_bai, mom_bam, mom_bai]
            .branch{ meta, kid_bam, kid_bai, dad_bam, dad_bai, mom_bam, mom_bai ->
                is_trio: (dad_bam != null && dad_bai != null) && (mom_bam != null && mom_bai != null)
                    return tuple ( meta, kid_bam, kid_bai, dad_bam, dad_bai, mom_bam, mom_bai )
                no_trio: (dad_bam == null && dad_bai == null) || (mom_bam == null && mom_bai == null)
                    return tuple ( meta, kid_bam, kid_bai, [], [], [], [] )
            }
        .set{ch_samples}

        trio_kids = ch_samples.is_trio.map{[it[0], it[1], it[2]]}
        trio_dads = ch_samples.is_trio.map{[it[0], it[3], it[4]]}
        trio_moms = ch_samples.is_trio.map{[it[0], it[5], it[6]]}

        //non_trio_kids = ch_samples.no_trio.map{[it[0], it[1], it[2]]} // These can be someone elses mom or dad
        //non_trio_dads = ch_samples.no_trio.map{[it[0], it[3], it[4]]} // These should all be empty
        //non_trio_moms = ch_samples.no_trio.map{[it[0], it[5], it[6]]} // These should all be empty

        // Get 3 vcf files back...we should collect all
        // Deal with multiple Parent VCFs in nextflow/GLNexus by naming them child.child/paternal/maternal
        // and the person running multiple trios per family will have to deal with it later (me)
        // Do we even want to merge though?

        // Maaaybe do an is_parent check and if not, run DEEPVARIANT

    } else {
        trio_kids = Channel.empty()
        trio_dads = Channel.empty()
        trio_moms = Channel.empty()
    }

    // Does splitting BAMs and copying to node make sense to reduce IO?

    // Only one of these is run depending on params.variant_caller (when clause condition is defined in the conf/modules.config)
    DEEPVARIANT               ( ch_bam_bai, ch_fasta, ch_fai )
    PEPPER_MARGIN_DEEPVARIANT ( ch_bam_bai, ch_fasta, ch_fai )
    DEEPTRIO                  ( trio_kids, ch_fasta, ch_fai, trio_dads, trio_moms)

    // Collect VCFs
    ch_snp_calls_vcf  = ch_snp_calls_vcf.mix(DEEPVARIANT.out.vcf)
    ch_snp_calls_vcf  = ch_snp_calls_vcf.mix(PEPPER_MARGIN_DEEPVARIANT.out.vcf)
    ch_snp_calls_vcf  = ch_snp_calls_vcf.mix(DEEPTRIO.out.vcf)

    // Collect GVCFs
    ch_snp_calls_gvcf = ch_snp_calls_gvcf.mix(DEEPVARIANT.out.gvcf)
    ch_snp_calls_gvcf = ch_snp_calls_gvcf.mix(PEPPER_MARGIN_DEEPVARIANT.out.gvcf)
    ch_snp_calls_gvcf = ch_snp_calls_gvcf.mix(DEEPTRIO.out.gvcf)

    // TODO: This only works with DeepVariant for now (remove PEPPER_MARGIN_DEEPVARIANT/Deeptrio?)
    ch_snp_calls_gvcf_region = DEEPVARIANT.out.gvcf_region

    // Group gVCFs from all samples per region for GLNexus
    /*ch_snp_calls_gvcf_region
        .map{
            meta, gvcf, region ->
            [['id':region.replaceAll(/[:-]/, "_")], gvcf] // meta.id = region
            }
        .groupTuple()           // size = number of samples (unknown)
        .set{ glnexus_regions } // channel: [ val(meta.id), path(one_gvcf_per_sample) ]*/

    // GLNexus has to bulk load all gVCFS, providing a BED file does not help in this step
    // So to make things a bit quicker, if we have a BED file, then we can cut out those regions...

       TABIX_EXTRA_GVCFS(ch_extra_gvcfs)

    /*ch_extra_gvcfs
        .join(TABIX_EXTRA_GVCFS.out.tbi)
        .groupTuple()
        .combine(SPLIT_BED_CHUNKS.out.split_beds)
        .view()*/

    BCFTOOLS_VIEW_REGIONS(ch_extra_gvcfs.join(TABIX_EXTRA_GVCFS.out.tbi).groupTuple(), ch_bed.collect())

    // This will add a "full" extra gVCF to all regions if not using a BED..
    // So maybe we should just concat the DV gVCFs together for each sample,
    // Then run GLNexus once...
    // If merging an already merged file with DV-files is faster, we can merge extra gVCFs first one time

    TABIX_DV(ch_snp_calls_gvcf)

    BCFTOOLS_CONCAT_DV(ch_snp_calls_gvcf.groupTuple().join(TABIX_DV.out.tbi.groupTuple()))

    BCFTOOLS_SORT_DV(BCFTOOLS_CONCAT_DV.out.vcf)
    // Concat DV and extra gvCFs together -> send to glnexus
    BCFTOOLS_SORT_DV.out.vcf
        .concat(BCFTOOLS_VIEW_REGIONS.out.vcf)
        .map { meta, gvcf -> [ ['id':'multisample'], gvcf ]}
        .groupTuple()
        .set{ ch_glnexus_in }

    ch_glnexus_in.view()

    // Then run GlNexus to join-call genotypes together
    GLNEXUS( ch_glnexus_in, ch_bed.collect() )

    // Tabix multisample bcf
    TABIX_CONCAT_UNMERGED_MULTISAMPLE(GLNEXUS.out.bcf)

    // Get versions
    ch_versions     = ch_versions.mix(DEEPVARIANT.out.versions)
    ch_versions     = ch_versions.mix(PEPPER_MARGIN_DEEPVARIANT.out.versions)
    ch_versions     = ch_versions.mix(DEEPTRIO.out.versions)
    ch_versions     = ch_versions.mix(GLNEXUS.out.versions)
    //ch_versions     = ch_versions.mix(BCFTOOLS_VIEW_REGIONS.out.versions)
    //ch_versions     = ch_versions.mix(TABIX_EXTRA_GVCFS.out.versions)
    //ch_versions     = ch_versions.mix(BCFTOOLS_CONCAT_SINGLESAMPLE.out.versions)
    //ch_versions     = ch_versions.mix(BCFTOOLS_SORT_CONCAT_SINGLESAMPLE.out.versions)
    //ch_versions     = ch_versions.mix(TABIX_CONCAT_UNMERGED_MULTISAMPLE.out.versions)
    //ch_versions     = ch_versions.mix(TABIX_CONCAT_UNMERGED_SINGLESAMPLE.out.versions)
    //ch_versions     = ch_versions.mix(TABIX_CONCAT_MERGED_SINGLESAMPLE.out.versions)


    emit:
    snp_calls_vcf  = BCFTOOLS_SORT_DV.out.vcf
    combined_bcf   = GLNEXUS.out.bcf
    versions       = ch_versions
}
