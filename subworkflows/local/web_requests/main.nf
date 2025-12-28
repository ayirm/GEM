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

workflow WEB_REQUESTS {
  take:
    gbk_ch // tuple val(meta), path(gbk_json)

  main:
    mapping_script = params.mapping_script
    mapping_json = UNIPROT_MAPPING(gbk_ch, mapping_script) // tuple val(meta), path(mapping_json) also version.yml

    goTerms_script = params.goTerms_script
    go_json = GO_TERM_FINDER(gbk_ch, goTerms_script)

  emit:
    enriched_json
}