process FIND_ANNOTATED_FAMILIES_BY_STRUCTURE {
    publishDir "${params.outDir}/foldseek", mode: "copy"
    label "general"

    input:
    path(m8s)

    output:
    path("annotated.txt"), optional: true, emit: annotated

    script:
    """
    awk '{print \$1}' ${m8s} | cut -d'-' -f 1 | sort | uniq > annotated.txt
    """
}
