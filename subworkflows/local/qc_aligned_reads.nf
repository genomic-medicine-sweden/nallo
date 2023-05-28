include { CRAMINO  } from '../../modules/local/cramino'
include { MOSDEPTH } from '../../modules/nf-core/mosdepth/main.nf'

workflow QC_ALIGNED_READS {

    take:
    ch_bam   // channel: [ val(meta), bam ]
    ch_bai   // channel: [ val(meta), bai ]
    ch_fasta // channel: [ val(meta), fasta ]

    main:
    ch_versions = Channel.empty()

    // Prepare inputs
    ch_bam
        .join(ch_bai)
        .set{ ch_cramino_in }
    
    // Prepare inputs
    ch_cramino_in
        .map{ meta, bam, bai -> [ meta, bam, bai, [] ] }
        .set{ ch_mosdepth_in }   
    
    CRAMINO (ch_cramino_in)
    MOSDEPTH(ch_mosdepth_in, ch_fasta)
    
    // Gather versions
    ch_versions = ch_versions.mix(CRAMINO.out.versions.first())
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}

