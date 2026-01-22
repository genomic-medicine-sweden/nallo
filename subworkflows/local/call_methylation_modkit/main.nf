include { MODKIT_PILEUP            } from '../../../modules/nf-core/modkit/pileup/main'
include { MODKIT_BEDMETHYLTOBIGWIG } from '../../../modules/nf-core/modkit/bedmethyltobigwig/main'
include { TABIX_TABIX              }  from '../../../modules/nf-core/tabix/tabix/main'
workflow CALL_METHYLATION_MODKIT {

    take:
    ch_bam_bai // channel: [ val(meta), bam, bai ]
    ch_fasta   // channel: [ val(meta), fasta ]
    ch_fai     // channel: [ val(meta), fai ]
    ch_bed     // channel: [ val(meta), bed ]
    modcodes   // String or List

    main:
    // Performs pileups per haplotype if the phasing workflow is on, set in config
    MODKIT_PILEUP(
        ch_bam_bai,
        ch_fasta.combine(ch_fai.map { _meta, fai -> fai }),
        ch_bed,
    )

    MODKIT_PILEUP.out.bedgz
        .transpose()
        .tap { ch_bedmethyl }
        // Only convert files with content
        .filter { _meta, bed -> gzNotEmptyBySize(bed) }
        .set { ch_bedmethyl_to_bigwig_in }

    TABIX_TABIX(
        ch_bedmethyl,
    )

    MODKIT_BEDMETHYLTOBIGWIG(
        ch_bedmethyl_to_bigwig_in,
        ch_fai,
        modcodes
    )

    emit:
    bed      = ch_bedmethyl                    // channel: [ val(meta), path(bed) ]
    tbi      = TABIX_TABIX.out.index           // channel: [ val(meta), path(tbi) ]
    bigwig   = MODKIT_BEDMETHYLTOBIGWIG.out.bw // channel: [ val(meta), path(bw) ]
}

def gzNotEmptyBySize(file_path) {
    File gzipFile = file_path.toFile()

    // A valid gzip file is at least ~18 bytes (header + footer)
    if (gzipFile.length() < 18) {
        return false
    }

    gzipFile.withInputStream { inputStream ->
        // ISIZE is stored in the last 4 bytes of the gzip file
        long footerOffset = gzipFile.length() - 4
        inputStream.skip(footerOffset)

        def footerStream = new DataInputStream(inputStream)
        int uncompressedSize = footerStream.readInt()

        return uncompressedSize > 0
    }
}
