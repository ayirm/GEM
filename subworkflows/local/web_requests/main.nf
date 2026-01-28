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
    mapping_script = params.mapping_script
    mapping_json = UNIPROT_MAPPING(gbk_ch, mapping_script) // tuple val(meta), path(mapping_json) also version.yml

    goTerms_script = params.goTerms_script
    go_json = GO_TERM_FINDER(gbk_ch, goTerms_script)

    kegg_requests_script = params.kegg_requests_script
    kegg_json = KEGG_REQUESTS(mapping_json, kegg_requests_script)

    json_merging_script = params.json_merging_script
    // TODO: This could result in an error since all the process channels also return a meta id, just double check it
    merged_json = JSON_MERGING(gbk_ch, mapping_json, go_json, kegg_json, json_merging_script) 

  emit:
    merged_json  // This goes to modelling workflow
}