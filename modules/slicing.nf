process SLICE {
    publishDir "${params.outDir}slices", mode: "copy"

    input:
    path fasta_path
    path annotations_path
    val min_slice_length

    output:
    path "sliced_fasta.csv"

    script:
    """
    python3 ${params.scriptDir}slice_proteins.py ${fasta_path} ${annotations_path} ${min_slice_length} sliced_fasta.csv
    """
}
