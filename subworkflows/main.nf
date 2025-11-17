include { GENS_OUTPUT } from './gens_output'

params.outdir = params.outdir ?: 'results'

workflow MAIN {

    def bam_fp = "/home/jakob/data/hg002_chr21_subsamp01_subset.bam"
    def bam_bai_fp = "/home/jakob/data/hg002_chr21_subsamp01_subset.bam.bai"
    def gvcf_fp = "/home/jakob/data/hg002_chr21.dnascope.gvcf.gz"
    def gvcf_tbi_fp = "/home/jakob/data/hg002_chr21.dnascope.gvcf.gz.tbi"
    def baf_positions_fp = ""

    def meta = [ id: 'SAMPLE1', sex: 'M', cohort: 'test' ]

    channel.of(
        tuple(
            meta,
            file(bam_fp),
            file(bam_bai_fp)
        )
    ).set { ch_bam }

    channel.of(
        tuple(
            meta,
            file(gvcf_fp),
            file(gvcf_tbi_fp),
            file(baf_positions_fp),
        )
    ).set { ch_gvcf }

    GENS_OUTPUT(
        ch_bam,
        ch_gvcf
    )

    ch_cov = GENS_OUTPUT.out.cov_bed_tbi
    ch_baf = GENS_OUTPUT.out.baf_bed_tbi

    ch_cov_baf = ch_cov.join(ch_baf)

    PUBLISH_GENS_OUTPUT(ch_cov_baf)

    ch_cov.view { it -> "cov_bed_tbi: $it" }
    ch_baf.view { it -> "baf_bed_tbi: $it" }
    GENS_OUTPUT.out.versions.view { "versions: $it" }
}

process PUBLISH_GENS_OUTPUT {
    publishDir "${params.outdir}/gens", mode: "copy"

    input:
    tuple val(meta), path(cov_bed_tbi), path(baf_bed_tbi)

    script:
    """
    ln -s $cov_bed_tbi
    ln -s $baf_bed_tbi
    """
}
