process EGGNOG_MAPPER {
    publishDir "${params.outDir}annotations/eggnog", mode: "copy"
    label "eggnog"

    input:
    path fasta_file
    path eggnog_data_dir
    path eggnog_diamond_db
    path eggnog_db

    output:
    path "${fasta_file.baseName}.emapper.hits"

    script:
    """
    emapper.py \
    --cpu ${task.cpus} \
    -i ${fasta_file} \
    --data_dir ${eggnog_data_dir} \
    -m diamond \
    --dmnd_db ${eggnog_diamond_db} \
    --database ${eggnog_db} \
    --output ${fasta_file.baseName}
    """
}