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
    publishDir "${params.outdir}/input", mode: "copy"
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

process EXTRACT_FIRST_STOCKHOLM_SEQUENCES {
    label "general"

    input:
    path msa
    path ids
    val mode

    output:
    path "family_reps.fasta"

    script:
    """
    python3 ${params.scriptDir}/extract_first_stockholm_sequences.py ${msa} ${ids} ${mode} family_reps.fasta
    """
}

process EXTRACT_UNIQUE_IDS {
    label "general"

    input:
    path tsvFiles

    output:
    file "unique_ids.txt"

    script:
    """
    for file in ${tsvFiles}; do 
        awk -F '\t' '{ print \$1 }' \$file | sort | uniq; 
    done | sort | uniq > unique_ids.txt
    """
}

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
    label "general"
    
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
    if [[ ! -s annotated_ids.txt ]]; then
        # If annotated_ids.txt is empty or does not exist, copy all lines from rep_names.txt to unannotated_ids.txt
        cp ${rep_names_txt} unannotated_ids.txt
    else
        awk 'NR==FNR {exclude[\$0]=1; next} !exclude[\$0]' annotated_ids.txt ${rep_names_txt} > unannotated_ids.txt
    fi
    """
}

process FIND_FASTA_BY_ID {
    label "general"

    input:
    path ids
    path fasta

    output:
    file "wp1_unannotated.fasta"

    script:
    """
    python3 ${params.scriptDir}/find_fasta_by_id.py ${ids} ${fasta} wp1_unannotated.fasta
    """
}
