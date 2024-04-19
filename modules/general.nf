process UNZIP_BZ2 {
    label "general"

    input:
    path bz2

    output:
    file "${bz2.baseName}"

    script:
    """
    bzip2 -d < ${bz2} > ${bz2.baseName}
    """
}

process REMOVE_HEADER {
    // publishDir "${params.outdir}/input", mode: "copy" // better link next subworkflow to work folder instead
    label "general"

    input:
    path file

    output:
    path "${file.baseName}_no_header.${file.extension}"

    script:
    """
    tail -n +2 ${file} > ${file.baseName}_no_header.${file.extension}
    """
}
