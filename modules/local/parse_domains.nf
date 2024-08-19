process PARSE_DOMAINS {
    input:
    path query_results
    path refined_families

    output:
    path "domain_results"

    """
    parse_domains.py ${query_results} ${refined_families}
    """
}
