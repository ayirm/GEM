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