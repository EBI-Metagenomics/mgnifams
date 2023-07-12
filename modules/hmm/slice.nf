process SLICE {
    publishDir 'data/output/hmm/slices', mode: 'copy', saveAs: { filename ->
            def newFilename = filename.replaceAll("family.fa_mafft.fa.hmm_scan.tblout_", "")
            "${newFilename}"
        }

    input:
    path tblout

    output:
    path "${tblout}_sliced.csv"

    script:
    """
    ${baseDir}/bin/slice_tblout.py ${tblout} ${tblout}_sliced.csv
    """
}

process SCANSLICE {
    input:
    path slice
    path msa

    script:
    """
    ${baseDir}/bin/scan_slice.py ${slice} ${msa}
    """
}