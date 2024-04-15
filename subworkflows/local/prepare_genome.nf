include { SAMTOOLS_FAIDX           } from '../../modules/nf-core/samtools/faidx/main'
include { MINIMAP2_INDEX           } from '../../modules/nf-core/minimap2/index/main'
include { GUNZIP as GUNZIP_FASTA   } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_VEP_CACHE } from '../../modules/nf-core/untar/main'

workflow PREPARE_GENOME {

    take:
    fasta_in // channel: [ val(meta), fasta ]
    ch_vep_cache       // channel: [mandatory for annotation] [ path(cache) ]


    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.empty()

    fasta_file = fasta_in.map{meta, file -> file}

    // Will not catch cases where fasta is bgzipped
    if ( params.fasta.endsWith('.gz') ) {
        GUNZIP_FASTA(fasta_in)
            .gunzip
            .collect()
            .set{ch_fasta}

        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions.first())
    } else {
        fasta_in
            .set{ch_fasta}
    }

    SAMTOOLS_FAIDX ( ch_fasta )
    MINIMAP2_INDEX ( ch_fasta )

    UNTAR_VEP_CACHE (ch_vep_cache)

    UNTAR_VEP_CACHE.out.untar
        .map{meta, files -> [files]}
        .collect()
        .set { untarred_vep }

    // Gather versions
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())
    ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

    emit:
    fasta         = ch_fasta                           // channel: [ val(meta), fasta ]
    fai           = SAMTOOLS_FAIDX.out.fai.collect()   // channel: [ val(meta), fai ]
    vep_resources = untarred_vep                       // channel: [ path(cache) ]
    mmi           = MINIMAP2_INDEX.out.index.collect() // channel: [ val(meta), mmi ]
    versions      = ch_versions                        // channel: [ versions.yml ]
}
