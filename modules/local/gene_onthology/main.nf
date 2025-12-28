process GO_TERM_FINDER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'python:3.11.9' }"

    input:
    tuple val(meta), path(parsed_gbk)
    path(python_script)

    output:
    tuple val(meta), path("*.json"), emit: go_terms
    path "versions.yml", emit: versions

    script:
    """
    python3 ${python_script} \
        --input_json ${parsed_gbk} \
        --go_json go_terms.json
    """
}