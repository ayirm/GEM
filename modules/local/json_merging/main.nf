process JSON_MERGING {
    tag "${meta.id}"
    label 'process_low'
    // irem bunu lowest yap sen

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/library/python:3.11.9' }"

    input:
    tuple val(meta), path(parsed_json)
    path(mapping_json)
    path(go_json)
    path(kegg_json)
    path (python_script)

    output:
    tuple val(meta), path("*_merged.json"), emit:merged
    path "versions.yml", emit: versions

    script:
    """
    python3 ${python_script} \
        --parsed_json ${parsed_json} \
        --mapping_json ${mapping_json} \
        --go_json ${go_json} \
        --kegg_json ${kegg_json} \
        --merged_json ${meta.id}_merged.json
    """
}