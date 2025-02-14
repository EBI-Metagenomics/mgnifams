process FIND_ANNOTATED_FAMILIES_BY_STRUCTURE {
    tag "$meta.id"
    label "general"

    input:
    tuple val(meta), path(m8s)

    output:
    tuple val(meta), path("annotated_structures.txt"), optional: true, emit: annotated_structures

    script:
    """
    awk '{print \$1}' ${m8s} | cut -d'-' -f 1 | sort | uniq > annotated_structures.txt
    """
}
