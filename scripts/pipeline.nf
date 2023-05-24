#!/usr/bin/env nextflow

// Define base path
params.base_path = "../../../"

// Define Python script path
params.db_to_fasta_script = "${params.base_path}db_to_fasta.py"

// A process for fetching data from a PostgreSQL database and writing it to a FASTA file
process fetchData {
    output:
    file("pgsql.fasta")

    shell:
    """
    python ${params.db_to_fasta_script}
    """
}

workflow {
    fetchData()
}
