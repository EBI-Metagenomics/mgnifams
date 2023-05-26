#!/usr/bin/env nextflow

// Define base path
params.base_path = "../../../"

// Define Python script paths

params.base_url = "http://localhost:8000/sequence_explorer"

// Fetch data from PostgreSQL db
// and write it to a FASTA file
// TODO change this with just input channel of a fasta file with the whole DB
process fetchData {
    output:
    path "pgsql.fasta"

    script:
    """
    db_to_fasta.py
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
    fasta_to_urls.py $pgsql ${params.base_url}
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
    query_urls.py $urls
    """
}

workflow {
    fetchData | createUrls | sendApiRequests
}