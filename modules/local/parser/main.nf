process GB_PARSER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/gb-parser:3.11' }"

    input:
    tuple val(meta), path(gb_file)
    path(python_script)

    output:
    tuple val(meta), path("*_parsed_gb.json"), emit: parsing
    path "versions.yml", emit: versions

    script:
    """
    python3 ${python_script} \
        --gb_file ${gb_file} \
        --parser_json ${meta.id}_parsed_gb.json
    """
}