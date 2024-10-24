include { GUNZIP as GUNZIP_FASTA   } from '../../modules/nf-core/gunzip/main'
include { MINIMAP2_INDEX           } from '../../modules/nf-core/minimap2/index/main'
include { SAMTOOLS_FAIDX           } from '../../modules/nf-core/samtools/faidx/main'
include { UNTAR as UNTAR_VEP_CACHE } from '../../modules/nf-core/untar/main'

workflow PREPARE_GENOME {

    take:
    fasta_in                   // channel: [mandatory] [ val(meta), path(fasta) ]
    gunzip_fasta               //    bool: should we gunzip fasta
    ch_vep_cache               // channel: [optional] [ val(meta), path(cache) ]
    split_vep_files            //    bool: are there vep extra files
    ch_vep_extra_files_unsplit // channel: [optional] [ val(meta), path(csv) ]

    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.empty()

    fasta_file = fasta_in.map{meta, file -> file}

    // Will not catch cases where fasta is bgzipped
    if ( gunzip_fasta ) {
        GUNZIP_FASTA ( fasta_in )
            .gunzip
            .collect()
            .set { ch_fasta }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions.first())
    } else {
        fasta_in
            .set { ch_fasta }
    }

    SAMTOOLS_FAIDX ( ch_fasta, [[],[]] )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    MINIMAP2_INDEX ( ch_fasta )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    UNTAR_VEP_CACHE (ch_vep_cache)
    ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

    UNTAR_VEP_CACHE.out.untar
        .map { meta, files -> [ files ] }
        .collect()
        .set { untarred_vep }

    // Read and store paths in the vep_plugin_files file
    if ( split_vep_files ) {
        ch_vep_extra_files_unsplit
            .splitCsv ( header:true )
            .map { row ->
                path = file(row.vep_files[0])
                if(path.exists()) {
                    return [path]
                } else {
                    error("\nVep database file ${path} does not exist.")
                }
            }
            .collect()
            .set { ch_vep_extra_files }
    } else {
        ch_vep_extra_files = Channel.value([])
    }

    emit:
    mmi             = MINIMAP2_INDEX.out.index.collect() // channel: [ val(meta), path(mmi) ]
    fai             = SAMTOOLS_FAIDX.out.fai.collect()   // channel: [ val(meta), path(fai) ]
    fasta           = ch_fasta                           // channel: [ val(meta), path(fasta) ]
    vep_resources   = untarred_vep                       // channel: [ path(cache) ]
    vep_extra_files = ch_vep_extra_files                 // channel: [ path(files) ]
    versions        = ch_versions                        // channel: [ versions.yml ]
}
