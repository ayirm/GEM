process GO_TERM_FINDER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11.9' :
        'docker.io/dogay/go-term-finder:3.11' }"

    input:
    tuple val(meta), path(parsed_gbk)
    path(python_script)

    output:
    tuple val(meta), path("*_go_terms.json"), emit: go_terms
    path "versions.yml", emit: versions

    script:
    """
    export XDG_CONFIG_HOME=\$PWD/.config
    export XDG_CACHE_HOME=\$PWD/.cache
    export MPLCONFIGDIR=\$PWD/.config/matplotlib
    mkdir -p \$XDG_CONFIG_HOME \$XDG_CACHE_HOME \$MPLCONFIGDIR
    python3 ${python_script} \
        --input_json ${parsed_gbk} \
        --go_json ${meta.id}_go_terms.json
    """
}