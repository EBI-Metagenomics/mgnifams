process DECOMPRESS_GZ {
    label "general"

    input:
    path gz_file

    output:
    file "${gz_file.baseName}"

    script:
    """
    pigz -dfk -p${task.cpus} ${gz_file} > ${gz_file.baseName}
    """
}
