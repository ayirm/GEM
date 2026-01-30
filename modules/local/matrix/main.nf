process GENERATE_MATRIX {
    tag "Matrix Generation"
    label 'process_lowest'
  
    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/library/python:3.11.9' }"

    input:
    path (json_input) // json file that contains the following: Pathway id, reaction id
    path (python_script)

    output:
    path "stoichiometric_matrix.xlsx", emit: matrix_xlsx

    script:
    """
    python3 ${python_script} \\
        --input ${json_input} \\
        --output stoichiometric_matrix.xlsx
    """
}
