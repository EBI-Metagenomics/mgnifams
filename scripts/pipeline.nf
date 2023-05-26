#!/usr/bin/env nextflow

// Define base path
params.base_path = "../../../"

// Define Python script paths
params.db_to_fasta_script = "${params.base_path}db_to_fasta.py"
params.fasta_to_urls_script = "${params.base_path}fasta_to_urls.py"
params.query_urls_script = "${params.base_path}query_urls.py"

params.base_url = "http://localhost:8000/sequence_explorer"

// Fetch data from PostgreSQL db
// and write it to a FASTA file
// TODO change this with just input channel of a fasta file with the whole DB
process fetchData {
    output:
    path "pgsql.fasta"

    script:
    """
    python3 ${params.db_to_fasta_script}
    """
}

// From protein ids, create urls for Django API
process createUrls {
    input:
    path pgsql

    output:
    path "urls.txt"

    script:
    """
    python3 ${params.fasta_to_urls_script} $pgsql ${params.base_url}
    """
}

// Hit Django API for every protein and save the response HTML page
process sendApiRequests {
    input:
    path urls

    output:
    path "*_page.html"

    script:
    """
    python3 ${params.query_urls_script} $urls
    """
}

workflow {
    fetchData | createUrls | sendApiRequests
}