process SLICE {
    publishDir 'data/output/slices', mode: 'copy', saveAs: { filename ->
            def newFilename = filename.replaceAll("_family_domains.tblout", "")
            "${newFilename}"
        }

    input:
    path tblout

    output:
    path "${tblout}_sliced.csv"

    script:
    """
    ${baseDir}/bin/slice_tblout.py ${tblout} ${tblout}_sliced.csv > /dev/null 2>&1
    """
}

process SCANSLICE {
    input:
    path slice
    path msa

    script:
    """
    ${baseDir}/bin/scan_slice.py ${slice} ${msa} > /dev/null 2>&1
    """
}