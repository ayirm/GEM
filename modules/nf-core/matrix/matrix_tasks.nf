process GENERATE_MATRIX {
    tag "Matrix Generation"
    cpus 1
    memory '1 GB'  // 12 GB yerine 1 GB yeter de artar bile
    time '10m'
  
    // Çıktıyı sonuçlar klasörüne kopyalar
    publishDir "${params.outdir}/metabolic_matrix", mode: 'copy'

    input:
    path json_input // Arkadaşının KEGG modülünden gelecek olan JSON

    output:
    path "stoichiometric_matrix.xlsx", emit: matrix_xlsx

    script:
    """
    python3 ${projectDir}/scripts/matrix/matrix_gen.py \\
        --input ${json_input} \\
        --output stoichiometric_matrix.xlsx
    """
}
