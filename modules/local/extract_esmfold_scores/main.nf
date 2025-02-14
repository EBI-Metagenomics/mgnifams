process EXTRACT_ESMFOLD_SCORES {
    tag "$meta.id"
    label "venv"

    input:
    tuple val(meta), path(scores)

    output:
    tuple val(meta), path("${meta.id}_scores.csv")                 , emit: csv
    tuple val(meta), path("${meta.id}_high_quality_structures.txt"), emit: txt

    script:
    """
    csv_file="${meta.id}_scores.csv"
    #name,length,plddt,ptm
    grep -E 'pLDDT' *_scores.txt | awk -F' |, ' '{print \$11 "," \$14 "," \$16 "," \$18}' >> "\$csv_file"

    txt_file="${meta.id}_high_quality_structures.txt"
    awk -F',' '\$3 >= 70 {print \$1}' ${meta.id}_scores.csv > "\$txt_file"
    """
}
