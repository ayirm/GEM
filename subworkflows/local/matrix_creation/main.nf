/*
Needs to get the following values:
- KEGG_REACTIONS (for RXN_ID)
- KEGG_PATHWAYS (for Compound ID)
then create a S matrix for return

It needs to have these
rxn_id, reaction string with(compund ids), reversibilty(as an extra guardrail), pathway info(for subsystem info)
*/

include { GENERATE_MATRIX } from "../../../modules/local/matrix/"

workflow MATRIX_GENERATION_WF {
    take:
    enriched_json // Girdi kanalı

    main:
    python_script = file(params.matrix_script)
    GENERATE_MATRIX(enriched_json, python_script)

    emit:
    matrix = GENERATE_MATRIX.out.matrix_xlsx
}
