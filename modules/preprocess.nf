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

process DECOMPRESS_BZ2 {
    label "general"

    input:
    path bz2_file

    output:
    file "${bz2_file.baseName}"

    script:
    """
    bzip2 -d < ${bz2_file} > ${bz2_file.baseName}
    """
}

process REMOVE_HEADER {
    label "general"

    input:
    path csv_file

    output:
    path "${csv_file.baseName}_no_header.${csv_file.extension}"

    script:
    """
    tail -n +2 ${csv_file} > ${csv_file.baseName}_no_header.${csv_file.extension}
    """
}
