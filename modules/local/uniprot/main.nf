process UNIPROT_MAPPING {
    tag "${meta.id}"
    label 'process_lowest'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/uniprot-mapping:3.11' }"

    input:
    tuple val(meta), path(parsed_json)
    path(python_script)

    output:
    tuple val(meta), path("*_mapping.json"), emit: mapping
    path "versions.yml", emit: versions

    script:
    """
    python3 ${python_script} \
        --parse_json ${parsed_json} \
        --mapping_json ${meta.id}_mapping.json
    """
}