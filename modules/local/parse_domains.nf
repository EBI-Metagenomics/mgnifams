process PARSE_DOMAINS {
    tag "$meta.id"
    
    input:
    tuple val(meta) , path(query_results)
    tuple val(meta2), path(pfam_mapping)
    tuple val(meta3), path(refined_families)

    output:
    tuple val(meta), path("domain_results")

    script:
    """
    parse_domains.py \\
        ${query_results} \\
        ${pfam_mapping} \\
        ${refined_families}
    """
}
