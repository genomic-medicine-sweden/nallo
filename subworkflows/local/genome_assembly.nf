include { HIFIASM                       } from '../../modules/nf-core/hifiasm'
include { YAK as YAK_PATERNAL           } from '../../modules/local/yak'
include { YAK as YAK_MATERNAL           } from '../../modules/local/yak'
include { GFASTATS as GFASTATS_MATERNAL } from '../../modules/nf-core/gfastats/main'
include { GFASTATS as GFASTATS_PATERNAL } from '../../modules/nf-core/gfastats/main'

workflow ASSEMBLY {

    take:
    ch_reads // channel: [ val(meta), fastq ]
    ch_ped   // channel: [ val(metameta) ]
    
    main:
    ch_versions = Channel.empty()

    if(params.hifiasm_mode == 'hifi-only') {
        HIFIASM ( ch_reads, [[],[]], [[],[]], [], []) // [ [meta], fastq ]
    } else if(params.hifiasm_mode == 'trio-binning') {
        ch_ped
            .combine(ch_reads
                .map{ meta, reads -> [meta.id, reads]
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
            .map{[it[0], it[2], it[4], it[6]]} // [meta, kid_reads, dad_reads, mom_reads]
            .branch{ meta, kid, dad, mom -> 
                is_trio: dad != null && mom != null
                    return tuple ( meta, kid, dad, mom ) 
                no_trio: dad == null || mom == null
                    return tuple ( meta, kid, [], [] )
            }
        .set{ch_samples}

        trio_kids = ch_samples.is_trio.map{[it[0], it[1]]}
        trio_dads = ch_samples.is_trio.map{[it[0], it[2]]}
        trio_moms = ch_samples.is_trio.map{[it[0], it[3]]}

        YAK_PATERNAL( trio_dads )
        YAK_MATERNAL( trio_moms )
        
        non_trio_kids = ch_samples.no_trio.map{[it[0], it[1]]} // These can be someone elses mom or dad 
        non_trio_dads = ch_samples.no_trio.map{[it[0], it[2]]} // These should all be empty
        non_trio_moms = ch_samples.no_trio.map{[it[0], it[3]]} // These should all be empty
        
        // TODO: This assembles the non-trio samples as well, add flag to disable?
        paternal_yak_or_empty = YAK_PATERNAL.out.yak.concat(non_trio_dads)
        maternal_yak_or_empty = YAK_MATERNAL.out.yak.concat(non_trio_moms)
        
        all_kid_reads = trio_kids.concat(non_trio_kids)
        // TODO: How to be really sure dad/mom are inputed correctly? Since test dataset all have same parents this needs to be double checked..
        
        HIFIASM (all_kid_reads, paternal_yak_or_empty, maternal_yak_or_empty, [], [] ) 
        
        ch_versions = ch_versions.mix(YAK_PATERNAL.out.versions.first())
        ch_versions = ch_versions.mix(YAK_MATERNAL.out.versions.first())
        
    }
    // Not the cleanest way, but better than to rely on hap_* in file names..
    GFASTATS_PATERNAL( HIFIASM.out.paternal_contigs,'fasta', '', '', [], [], [], [] )
    GFASTATS_MATERNAL( HIFIASM.out.maternal_contigs,'fasta', '', '', [], [], [], [] )

    GFASTATS_PATERNAL.out.assembly.combine(GFASTATS_MATERNAL.out.assembly, by: 0).set{ ch_dual_assembly_fa }
    
    ch_versions = ch_versions.mix(HIFIASM.out.versions.first())
    ch_versions = ch_versions.mix(GFASTATS_PATERNAL.out.versions.first())
    ch_versions = ch_versions.mix(GFASTATS_MATERNAL.out.versions.first())
    
    emit:
    assembled_haplotypes = ch_dual_assembly_fa // channel: [ [meta], paternal_fa, maternal_fa ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

