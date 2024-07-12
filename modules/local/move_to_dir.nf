process MOVE_TO_DIR {
    label "venv"

    input:
    path(files)
    val(dir_name)

    output:
    path("${dir_name}")

    script:
    """
    mkdir -p ${dir_name}
    cp -r ${files} ${dir_name}
    """
}
