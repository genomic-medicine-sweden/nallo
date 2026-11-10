include { GAWK as GAWK_DMR_REGIONS    } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_PILEUP_REGIONS } from '../../../modules/nf-core/gawk/main'
include { MODKIT_BEDMETHYLTOBIGWIG    } from '../../../modules/nf-core/modkit/bedmethyltobigwig/main'
include { MODKIT_DMR                  } from '../../../modules/nf-core/modkit/dmr/main'
include { MODKIT_PILEUP               } from '../../../modules/nf-core/modkit/pileup/main'
include { TABIX_TABIX                 } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_METHYLATION_MODKIT {
    take:
    ch_bam_bai   // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta     // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai       // channel: [mandatory] [ val(meta), path(fai) ]
    ch_bed       // channel: [optional]  [ val(meta), path(bed) ]
    ch_regions   // channel: [optional]  [ val(meta), path(tsv) ]
    modcodes     // string: [mandatory] modification codes for bigwigs
    val_skip_dmr // boolean: [mandatory] skip the HP1 vs HP2 comparison over ch_regions

    main:
    ch_fasta_fai = ch_fasta
        .combine(ch_fai.map { _meta, fai -> fai })
        .collect()

    // Trim BED to 3 columns; modkit pileup --include-bed requires a 3-column BED
    ch_bed_branches = ch_bed.branch { _meta, bed ->
        with_bed: bed
        no_bed: !bed
    }

    GAWK_PILEUP_REGIONS(
        ch_bed_branches.with_bed,
        [],
        false,
    )

    // .collect() converts from queue channel back to value channel so all samples get the bed
    ch_bed_for_modkit = ch_bed_branches.no_bed.mix(GAWK_PILEUP_REGIONS.out.output).collect()

    // Performs pileups per haplotype if the phasing workflow is on, set in config
    MODKIT_PILEUP(
        ch_bam_bai,
        ch_fasta_fai,
        ch_bed_for_modkit,
    )

    // Only convert files with content
    ch_bedmethyl_to_bigwig_in = MODKIT_PILEUP.out.bedgz
        .transpose()
        .tap { ch_bedmethyl }
        .filter { _meta, bed -> gzNotEmptyBySize(bed) }

    TABIX_TABIX(
        ch_bedmethyl
    )

    MODKIT_BEDMETHYLTOBIGWIG(
        ch_bedmethyl_to_bigwig_in,
        ch_fai,
        modcodes,
    )

    ch_dmr = channel.empty()

    //
    // Compare haplotype 1 with haplotype 2 over the methylation regions
    //
    if (!val_skip_dmr) {
        // One BED row per region and label; background models have one row per region and data category
        GAWK_DMR_REGIONS(
            ch_regions,
            [],
            false,
        )

        // Pair each bedMethyl with its index by file name, then pair hp1 with hp2 per sample
        ch_bedmethyl_tbi = ch_bedmethyl
            .map { meta, bed -> [bed.simpleName, meta, bed] }
            .join(TABIX_TABIX.out.index.map { _meta, tbi -> [tbi.simpleName, tbi] })
            .map { _name, meta, bed, tbi -> [meta, bed, tbi] }
            .branch { _meta, bed, _tbi ->
                hp1: bed.name.endsWith('_hp1.bed.gz')
                hp2: bed.name.endsWith('_hp2.bed.gz')
            }

        ch_dmr_input = ch_bedmethyl_tbi.hp1
            .join(ch_bedmethyl_tbi.hp2, failOnMismatch: true, failOnDuplicate: true)
            .multiMap { meta, hp1_bed, hp1_tbi, hp2_bed, hp2_tbi ->
                hp1: [meta, hp1_bed, hp1_tbi]
                hp2: [meta, hp2_bed, hp2_tbi]
            }

        MODKIT_DMR(
            ch_dmr_input.hp1,
            ch_dmr_input.hp2,
            GAWK_DMR_REGIONS.out.output.collect(),
            ch_fasta_fai,
        )

        ch_dmr = MODKIT_DMR.out.bed
    }

    emit:
    bed    = ch_bedmethyl                    // channel: [ val(meta), path(bed) ]
    tbi    = TABIX_TABIX.out.index           // channel: [ val(meta), path(tbi) ]
    bigwig = MODKIT_BEDMETHYLTOBIGWIG.out.bw // channel: [ val(meta), path(bw) ]
    dmr    = ch_dmr                          // channel: [ val(meta), path(bed) ]
}

def gzNotEmptyBySize(file_path) {
    def gzipFile: File = file_path.toFile()
    // When modkit produces an emty file, its size seems to be 168 bytes
    if (gzipFile.length() > 168) {
        return true
    }
}
