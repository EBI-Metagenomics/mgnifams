process INTERPROSCAN {
    label 'ips'

    container 'quay.io/microbiome-informatics/genomes-pipeline.ips:5.62-94.0'
    containerOptions '--bind data:/opt/interproscan-5.62-94.0/data'

    input:
    file faa_fasta
    path interproscan_db

    output:
    path '*.tsv', emit: tsv

    script:
    """
    interproscan.sh \\
        -cpu ${task.cpus} \\
        --input ${faa_fasta} \\
        -f TSV \\
        -dp \\
        -appl TIGRFAM,SUPERFAMILY,PANTHER,Gene3D,Hamap,ProSiteProfiles,CDD,AntiFam,Pfam \\
        -o ${faa_fasta.baseName}.tsv
    """
}
