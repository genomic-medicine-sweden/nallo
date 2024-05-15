include { SOMALIER_EXTRACT } from '../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE  } from '../../modules/nf-core/somalier/relate/main'

workflow BAM_ASSIGN_SEX {

    take:
    ch_bam_bai        // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta          // channel: [ val(meta), path(fasta) ]
    ch_fai            // channel: [ val(meta), path(fai) ]
    ch_somalier_sites // channel: [ val(meta), path(somalier_sites_vcf) ]
    ch_ped            // channel: [ val(meta), path(ped) ]

    main:

    ch_versions = Channel.empty()

    // Extract sites
    SOMALIER_EXTRACT (
        ch_bam_bai,
        ch_fasta.map { it[1] },
        ch_fai.map { it[1] },
        ch_somalier_sites.map { it[1] }
    )
    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions)

    SOMALIER_EXTRACT.out.extract
        .map { meta, extract -> [ [ id: 'multisample' ], extract ] }
        .groupTuple()
        .join( ch_ped.map { ped -> [ [ id:'multisample'], ped ] } )
        .set { ch_somalier_relate_in }

    // Infer sex
    SOMALIER_RELATE ( ch_somalier_relate_in, [] )
    ch_versions = ch_versions.mix(SOMALIER_RELATE.out.versions)

    SOMALIER_RELATE.out.samples_tsv
        .map { meta, tsv -> tsv }
        .splitCsv(header: true, sep: '\t')
        .set { somalier_tsv }

    // Hard error if sex could not be inferred for samples
    // with sex not set in samplesheet (sex == 0)
    somalier_tsv
        .map { it ->
            assert (it.sex.toInteger() != 0): error ("ERROR: Sex could not be automatically inferred for $it.sample_id. Please inspect manually and set sex samplesheet.")
        }

    somalier_tsv
        .map { it ->
            new_meta = [
                sample_id: it.sample_id,
                infered_sex: it.sex,
            ]
            [ it.sample_id, new_meta ]
        }
        .set { ch_somalier_sex }

    ch_bam_bai
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .join( ch_somalier_sex )
        .map { id, meta, bam, bai, meta2 ->
            new_meta = [
                id          : meta.id,
                family_id   : meta.family_id,
                paternal_id : meta.paternal_id,
                maternal_id : meta.maternal_id,
                sex         : meta2.infered_sex.toInteger(), // Use sex from somalier, only updated if sex is unknown
                phenotype   : meta.phenotype,
                single_end  : meta.single_end
            ]
            [ new_meta, bam, bai ]
        }
        .set { ch_updated_sex }

    emit:
    bam              = ch_updated_sex.map { meta, bam, bai -> [ meta, bam ] } // channel: [ val(meta), path(bam) ]
    bai              = ch_updated_sex.map { meta, bam, bai -> [ meta, bai ] } // channel: [ val(meta), path(bai) ]
    bam_bai          = ch_updated_sex                                         // channel: [ val(meta), path(bam), path(bai) ]
    somalier_samples = SOMALIER_RELATE.out.samples_tsv                        // channel: [ val(meta), path(samples_tsv) ]
    somalier_pairs   = SOMALIER_RELATE.out.pairs_tsv                          // channel: [ val(meta), path(pairs_tsv) ]
    versions = ch_versions                                                    // channel: [ versions.yml ]
}

