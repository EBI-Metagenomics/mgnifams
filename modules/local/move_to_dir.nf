process MOVE_TO_DIR {
    label "venv"

    input:
    tuple val(meta), path(files)
    val(dir_name)

    output:
    tuple val(meta), path("${dir_name}")

    script:
    """
    mkdir -p ${dir_name}
    cp -r ${files} ${dir_name}
    """
}
