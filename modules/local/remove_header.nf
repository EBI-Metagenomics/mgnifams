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
