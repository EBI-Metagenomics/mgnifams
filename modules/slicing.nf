process SLICE {
    publishDir 'data/output/slices', mode: 'copy', saveAs: { filename ->
        def newFilename = filename.replaceAll("_family_domains.tblout", "")
        "${newFilename}"
    }

    input:
    path tblout

    output:
    path "${tblout}_slices.csv"

    script:
    """
    ${baseDir}/bin/slice_tblout.py ${tblout} ${tblout}_slices.csv > /dev/null 2>&1
    """
}

process PRINT_SLICED {
    publishDir 'data/output/slices', mode: 'copy', saveAs: { filename ->
        def newFilename = filename.replaceAll("_msa.fa", "")
        "${newFilename}"
    }

    input:
    tuple path(msa), path(slice)

    output:
    path "${msa}_sliced.fa", optional: true

    script:
    """
    ${baseDir}/bin/print_sliced_seq_reps.py ${slice} ${msa} ${msa}_sliced.fa ${params.min_slice_length}
    """
}

process CONCAT_FASTA {
    publishDir 'data/output/slices', mode: 'copy'

    input:
    path fasta_files

    output:
    path "concatenated.fasta"

    """
    cat ${fasta_files.join(' ')} > concatenated.fasta
    """
}