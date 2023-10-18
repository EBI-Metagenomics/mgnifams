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

process FIND_UNANNOTATED_IDS {
    publishDir "${params.outdir}", mode: "copy"

    input:
    path annotations_csv
    path rep_names_txt

    output:
    path "unannotated_ids.txt"

    script:
    """
    # Extract the FamilyIDs from the CSV file and save them into a temporary file
    awk -F, 'NR > 1 {print \$1}' ${annotations_csv} | sort -u > annotated_ids.txt

    # Compare the two files and save the FamilyIDs from rep_names_txt that do not exist in annotated_ids.txt
    grep -vFf annotated_ids.txt ${rep_names_txt} > unannotated_ids.txt
    """
}
