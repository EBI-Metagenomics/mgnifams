process EXTRACT_ESMFOLD_SCORES {
    label "venv"
    publishDir "${params.outDir}/pdb", mode: "copy"

    input:
    tuple val(meta), path(scores)

    output:
    tuple val(meta), path("${meta.id}_scores.csv")                 , emit: csv
    tuple val(meta), path("${meta.id}_high_quality_structures.txt"), emit: txt

    script:
    """
    csv_file="${meta.id}_scores.csv"
    echo "name,length,pLDDT,pTM" > "\$csv_file"
    grep -E 'pLDDT' ${meta.id}_scores.txt | awk -F' |, ' '{print \$11 "," \$14 "," \$16 "," \$18}' >>  "\$csv_file"

    txt_file="${meta.id}_high_quality_structures.txt"
    awk -F',' 'NR>1 && \$3 >= 70 {print \$1}' ${meta.id}_scores.csv > "\$txt_file"
    """
}
