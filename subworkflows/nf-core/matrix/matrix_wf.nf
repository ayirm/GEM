include { GENERATE_MATRIX } from '../../../modules/nf-core/matrix/matrix_tasks.nf'

workflow MATRIX_GENERATION_WF {
    take:
    enriched_json // Girdi kanalı

    main:
    GENERATE_MATRIX(enriched_json)

    emit:
    matrix = GENERATE_MATRIX.out.matrix_xlsx
}
