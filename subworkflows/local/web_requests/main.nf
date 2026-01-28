/*
Needs to run the following
PARSE_GBK (either in this workflow or as an addition to the annotation workflow, needs to run a python script)
  └── gbk.json
        ├── UNIPROT_KEGG(id_mapping stage)
        │      └── mapping.json (kegg_stage)
        │             ├── KEGG_PATHWAYS
        │             ├── KEGG_REACTIONS
        │             └── KEGG_KO
        └── GO_MAPPING(go_terms stage)
               └── go.json

All of them needs to be merged into one giant .json file later, probably, i think sending one big file to parse is easier than sending multiple files to parse
*/

include { UNIPROT_MAPPING } from "../../../modules/local/uniprot"
include { GO_TERM_FINDER  } from "../../../modules/local/gene_onthology"
include { KEGG_REQUESTS   } from "../../../modules/local/kegg"
include { JSON_MERGING    } from "../../../modules/local/json_merging"

workflow WEB_REQUESTS {
  take:
    gbk_ch // tuple val(meta), path(gbk_json)

  main:
    mapping_script = channel.fromPath(params.mapping_script)
    UNIPROT_MAPPING(gbk_ch, mapping_script) // tuple val(meta), path(mapping_json) also version.yml

    goTerm_script = channel.fromPath(params.goTerm_script)
    GO_TERM_FINDER(gbk_ch, goTerm_script)

    kegg_requests_script = channel.fromPath(params.kegg_requests_script)
    KEGG_REQUESTS(UNIPROT_MAPPING.out.mapping, kegg_requests_script)

    json_merging_script = channel.fromPath(params.json_merging_script)
    JSON_MERGING(gbk_ch, UNIPROT_MAPPING.out.mapping, GO_TERM_FINDER.out.go_terms, KEGG_REQUESTS.out.kegg, json_merging_script) 

  emit:
    merged_json = JSON_MERGING.out.merged
    versions    = JSON_MERGING.out.versions
}