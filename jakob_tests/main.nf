include { PREPARE_GENS_INPUTS } from '../subworkflows/local/prepare_gens_inputs/main'

params.outdir = params.outdir ?: 'results'


workflow {

    def bam_fp = "/home/jakob/data/hg002_chr21_subsamp01_subset.bam"
    def bam_bai_fp = "/home/jakob/data/hg002_chr21_subsamp01_subset.bam.bai"
    def gvcf_fp = "/home/jakob/data/hg002_chr21.dnascope.gvcf.gz"
    def gvcf_tbi_fp = "/home/jakob/data/hg002_chr21.dnascope.gvcf.gz.tbi"
    def baf_positions_fp = "/home/jakob/data/gnomad_hg38.0.05.txt"
    def gatk_header_fp = "/home/jakob/data/header_tsv_gatk_mosdepth"
    def pon = "/fs2/viktor/LRS/gens/own_samples_mixed/50_100_hg38_5p_mosdepth.hdf5"

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

    PREPARE_GENS_INPUTS(
        ch_bam,
        ch_gvcf,
        gatk_header_fp,
        pon
    )

    ch_cov = PREPARE_GENS_INPUTS.out.cov_bed_tbi
    ch_baf = PREPARE_GENS_INPUTS.out.baf_bed_tbi

    ch_cov_baf = ch_cov.join(ch_baf)

    PUBLISH_GENS_INPUTS(ch_cov_baf)

    ch_cov.view { it -> "cov_bed_tbi: $it" }
    ch_baf.view { it -> "baf_bed_tbi: $it" }
    PREPARE_GENS_INPUTS.out.versions.view { it -> "versions: $it" }
}

process PUBLISH_GENS_INPUTS {
    publishDir "${params.outdir}/gens", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'ubuntu:22.04' }"

    input:
    tuple val(meta), path(cov_bed), path(cov_bed_tbi), path(baf_bed), path(baf_bed_tbi)

    script:
    """
    rm $cov_bed_tbi
    rm $baf_bed_tbi
    ln -s $cov_bed_tbi
    ln -s $baf_bed_tbi
    """

    stub:
    """
    touch $cov_bed_tbi
    touch $baf_bed_tbi
    """
}
