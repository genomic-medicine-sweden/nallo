include { GUNZIP as GUNZIP_FASTA   } from '../../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX           } from '../../../modules/nf-core/samtools/faidx/main'
include { UNTAR as UNTAR_VEP_CACHE } from '../../../modules/nf-core/untar/main'

workflow PREPARE_REFERENCES {
    take:
    fasta_in // channel: [ val(meta), path(fasta) ]
    fai_in // channel: [ val(meta), path(fai) ]
    ch_vep_cache // channel: [ val(meta), path(cache) ]
    gunzip_fasta // boolean: should we gunzip fasta
    untar_vep_cache // boolean: should we untar vep cache

    main:
    ch_fasta = channel.empty()

    // Will not catch cases where fasta is bgzipped
    if (gunzip_fasta) {
        ch_fasta = GUNZIP_FASTA(fasta_in).gunzip.collect()
    }
    else {
        ch_fasta = fasta_in
    }

    if (!fai_in) {
        SAMTOOLS_FAIDX(
            ch_fasta.map { meta, fasta -> [meta, fasta, []] },
            false,
        )

        ch_fai = SAMTOOLS_FAIDX.out.fai.collect()
    }
    else {
        ch_fai = fai_in
    }

    if (untar_vep_cache) {
        UNTAR_VEP_CACHE(
            ch_vep_cache
        )
    }

    emit:
    fai           = ch_fai // channel: [ val(meta), path(fai) ]
    fasta         = ch_fasta // channel: [ val(meta), path(fasta) ]
    vep_resources = untar_vep_cache ? UNTAR_VEP_CACHE.out.untar.collect() : ch_vep_cache // channel: [ val(meta), path(cache) ]
}
