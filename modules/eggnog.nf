process EGGNOG_MAPPER {
    publishDir "${params.outDir}annotations/eggnog", mode: "copy"
    container "quay.io/microbiome-informatics/genomes-pipeline.eggnog-mapper:v2.1.11"
    cpus 10

    input:
    path fasta_file
    path eggnog_diamond_db

    output:
    path "${fasta_file.baseName}.emapper.hits"

    script:
    """
    emapper.py \
    -i ${fasta_file} \
    --dmnd_db ${eggnog_diamond_db} \
    -m diamond \
    --cpu ${task.cpus} \
    --output ${fasta_file.baseName}
    """
    // python3 ${params.eggnog_dir}emapper.py -m diamond -i ${fasta_file} --cpu ${task.cpus} --output ${fasta_file.baseName} // old local, TODO remove when debugged
}