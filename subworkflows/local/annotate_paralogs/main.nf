include { PARAPHRASE } from '../../../modules/nf-core/paraphrase/main'

workflow ANNOTATE_PARALOGS {

    take:
    ch_json                  // channel: [ val(meta), path(json) ]
    paraphrase_output_format // string: Output format for paraphrase (json or tsv)
    ch_paraphrase_rules_yaml // channel: [ val(meta), path(yaml) ]

    main:
    /*
     * Paraphase emits one JSON per sample, but the JSON itself contains no sample identifier.
     * Therefore, sample names must be passed alongside the JSONs and grouped per family.
     *
     * The order of JSON files and sample names must remain aligned, since paraphrase assigns sample names to JSONs by positional order.
     */
    ch_json
        .map { meta, json -> [ [ 'id': meta.family_id ], json, meta.id ] }
        .groupTuple()
        .set { ch_paraphase_jsons_per_family }

    PARAPHRASE (
        ch_paraphase_jsons_per_family,
        ch_paraphrase_rules_yaml,
        paraphrase_output_format == 'tsv',
    )

    emit:
    json = PARAPHRASE.out.json // channel: [ val(meta), path(json) ]
    tsv  = PARAPHRASE.out.tsv // channel: [ val(meta), path(tsv) ]
}
