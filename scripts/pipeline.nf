#!/usr/bin/env nextflow

params.base_url = "http://localhost:8000/sequence_explorer"

Channel
    .fromPath("pgsql.fasta")
    .splitText(by: 50)
    .set { split_files }

process createUrls {
    input:
    val fasta_content

    output:
    path "urls.txt"

    script:
    """
    echo '${fasta_content}' | fasta_to_urls.py ${params.base_url}
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
    createUrls(split_files)
    sendApiRequests(createUrls.out)
}