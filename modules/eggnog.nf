process EGGNOG_MAPPER {
    publishDir 'data/output/annotations/eggnog', mode: 'copy'
    cpus 10

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}.emapper.hits"

    script:
    """
    python3 ${params.eggnog_dir}emapper.py -m diamond -i ${fasta_file} --cpu ${task.cpus} --output ${fasta_file.baseName}
    """
}