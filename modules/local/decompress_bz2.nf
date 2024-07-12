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
