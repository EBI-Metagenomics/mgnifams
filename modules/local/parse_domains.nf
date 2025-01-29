process PARSE_DOMAINS {
    input:
    tuple val(meta) , path(query_results)
    tuple val(meta2), path(refined_families)

    output:
    tuple val(meta), path("domain_results")

    script:
    """
    parse_domains.py ${query_results} ${refined_families}
    """
}
