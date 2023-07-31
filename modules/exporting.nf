process EXPORT_PROTEINS_CSV {
    publishDir 'data/output/tables', mode: 'copy'

    input:
    path fasta

    output:
    path "mgnifams_proteins.csv"

    """
    python3 ${baseDir}/bin/export_proteins_csv.py ${fasta} mgnifams_proteins.csv
    """
}

process EXPORT_CLUSTERING_CSV {
    publishDir 'data/output/tables', mode: 'copy'

    input:
    path cluster_file

    output:
    path "mgnifams_clusters.csv"

    """
    echo 'ProteinID,FamilyID' > mgnifams_clusters.csv && awk -v OFS="," '{print \$2,\$1}' ${cluster_file} >> mgnifams_clusters.csv
    """
}

process EXPORT_FAMILIES_CSV {
    publishDir 'data/output/tables', mode: 'copy'

    input:
    path rep_names
    val mode

    output:
    path "mgnifams_${mode}_families.csv"

    """
    python3 ${baseDir}/bin/export_families_csv.py ${rep_names} mgnifams_${mode}_families.csv ${mode}
    """
}

process EXPORT_KNOWN_ANNOTATIONS_CSV {
    input:
    path known_fasta

    output:
    path "${known_fasta.baseName}_annotations.csv", optional: true

    """
    python3 ${baseDir}/bin/export_known_annotations_csv.py ${known_fasta} ${known_fasta.baseName}_annotations.csv
    """
}

process CONCAT_ANNOTATIONS {
    publishDir 'data/output/tables', mode: 'copy'
    
    input:
    path known_fasta
    val mode

    output:
    path "mgnifams_${mode}_annotations.csv"

    """
    echo "FamilyID,Annotation,Description,Source,IsKnown" > mgnifams_${mode}_annotations.csv
    for csv_file in ${known_fasta}
    do
        cat \$csv_file >> mgnifams_${mode}_annotations.csv
    done
    """
}

process EXPORT_INTERPRO_ANNOTATIONS_CSV {
    input:
    path interpro

    output:
    path "interpro_annotations.csv", optional: true

    """
    python3 ${baseDir}/bin/export_interpro_annotations_csv.py ${interpro} interpro_annotations.csv
    """
}

process EXPORT_EGGNOG_ANNOTATIONS_CSV {
    input:
    path eggnog

    output:
    path "eggnog_annotations.csv", optional: true

    """
    python3 ${baseDir}/bin/export_eggnog_annotations_csv.py ${eggnog} eggnog_annotations.csv
    """
}

process EXPORT_UNIPROT_ANNOTATIONS_CSV {
    input:
    path uniprot

    output:
    path "${uniprot.baseName}_uniprot_annotations.csv", optional: true

    """
    python3 ${baseDir}/bin/export_uniprot_annotations_csv.py ${uniprot} ${uniprot.baseName}_uniprot_annotations.csv
    """
}