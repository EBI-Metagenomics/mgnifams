process CONCAT_FILES {
    input:
    path file1
    path file2

    output:
    path "combined_input.fa"

    """
    gunzip -c ${file1} > ${file1.baseName}
    cat ${file1.baseName} ${file2} > combined_input.fa
    """
}