process SLICE {
    publishDir "${params.outDir}slices", mode: "copy", saveAs: { filename ->
        def newFilename = filename.replaceAll("_family_domains.tblout", "")
        "${newFilename}"
    }

    input:
    path tblout

    output:
    path "${tblout}_slices.csv"

    script:
    """
    python3 ${params.scriptDir}slice_tblout.py ${tblout} ${tblout}_slices.csv > /dev/null 2>&1
    """
}

process PRINT_SLICED {
    publishDir "${params.outDir}slices", mode: "copy", saveAs: { filename ->
        def newFilename = filename.replaceAll("_msa.fa", "")
        "${newFilename}"
    }

    input:
    tuple path(msa), path(slice)

    output:
    path "${msa}_sliced.fa", optional: true

    script:
    """
    python3 ${params.scriptDir}print_sliced_seq_reps.py ${slice} ${msa} ${msa}_sliced.fa ${params.min_slice_length}
    """
}

process CONCAT_FASTA {
    publishDir "${params.outDir}slices", mode: "copy"

    input:
    path fasta_files

    output:
    path "concatenated.fasta"

    """
    cat ${fasta_files.join(' ')} > concatenated.fasta
    """
}