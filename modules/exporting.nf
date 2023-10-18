process EXPORT_PROTEINS_CSV {
    publishDir "${params.outdir}tables", mode: "copy"
    label "venv"

    input:
    path fasta

    output:
    path "mgnifams_proteins.csv"

    """
    python3 ${params.scriptDir}export_proteins_csv.py ${fasta} mgnifams_proteins.csv
    """
}

process EXPORT_CLUSTERING_CSV {
    publishDir "${params.outdir}tables", mode: "copy"

    input:
    path cluster_file

    output:
    path "mgnifams_clusters.csv"

    """
    echo 'ProteinID,FamilyID' > mgnifams_clusters.csv && awk -v OFS="," '{print \$2,\$1}' ${cluster_file} >> mgnifams_clusters.csv
    """
}

process EXPORT_FAMILIES_CSV {
    publishDir "${params.outdir}tables", mode: "copy"
    label "venv"
    
    input:
    path rep_names
    val mode

    output:
    path "mgnifams_${mode}_families.csv"

    """
    python3 ${params.scriptDir}export_families_csv.py ${rep_names} mgnifams_${mode}_families.csv ${mode}
    """
}

process EXPORT_KNOWN_ANNOTATIONS_CSV {
    label "venv"
    
    input:
    path known_fasta

    output:
    path "${known_fasta.baseName}_annotations.csv", optional: true

    """
    python3 ${params.scriptDir}export_known_annotations_csv.py ${known_fasta} ${known_fasta.baseName}_annotations.csv
    """
}

process EXPORT_BLASTP_ANNOTATIONS_CSV {
    label "venv"
    
    input:
    path fasta

    output:
    path "${fasta.baseName}_annotations.csv", optional: true

    """
    python3 ${params.scriptDir}export_blastp_annotations_csv.py ${fasta} ${fasta.baseName}_annotations.csv
    """
}

process CONCAT_ANNOTATIONS {
    publishDir "${params.outdir}tables", mode: "copy"
    
    input:
    path fasta
    val mode

    output:
    path "mgnifams_${mode}_annotations.csv"

    """
    echo "FamilyID,Annotation,Description,Source,IsKnown" > mgnifams_${mode}_annotations.csv
    for csv_file in ${fasta}
    do
        cat \$csv_file >> mgnifams_${mode}_annotations.csv
    done
    """
}

process EXPORT_INTERPRO_ANNOTATIONS_CSV {
    label "venv"
    
    input:
    path interpro

    output:
    path "interpro_annotations.csv", optional: true

    """
    python3 ${params.scriptDir}export_interpro_annotations_csv.py ${interpro} interpro_annotations.csv
    """
}

process EXPORT_EGGNOG_ANNOTATIONS_CSV {
    label "venv"

    input:
    path eggnog

    output:
    path "eggnog_annotations.csv", optional: true

    """
    python3 ${params.scriptDir}export_eggnog_annotations_csv.py ${eggnog} eggnog_annotations.csv
    """
}

process EXPORT_UNIPROT_ANNOTATIONS_CSV {
    label "venv"
    
    input:
    path uniprot

    output:
    path "${uniprot.baseName}_uniprot_annotations.csv", optional: true

    """
    python3 ${params.scriptDir}export_uniprot_annotations_csv.py ${uniprot} ${uniprot.baseName}_uniprot_annotations.csv
    """
}