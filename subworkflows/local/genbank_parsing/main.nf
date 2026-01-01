include { GB_PARSER } from "../../../modules/local/parser"

workflow GENBANK_FILE_PARSING {
    take:
      gb_file_ch // tuple val(meta), path(genbank_file)

    main:
      parser_script = params.parser_script
      gbk_json = GB_PARSER(gb_file_ch, parser_script)
    
    emit:
      gbk_json
}