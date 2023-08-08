process CONCAT_FILES {
    input:
    path mgy_folder
    path uniprot_sp

    output:
    path "combined_input.fa"

    script:
    """
    for file in ${mgy_folder}; do
        baseName=\${file%.gz}
        gunzip -c \$file > \$baseName
        cat \$baseName >> combined_input.fa
    done
    cat ${uniprot_sp} >> combined_input.fa
    """
}