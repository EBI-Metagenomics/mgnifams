process COMBINE_HH_RESULTS {
    label "general"

    input:
    tuple val(meta), path(hhr_folder)

    output:
    tuple val(meta), path("self_hits.tsv"), optional: true, emit: pfam_hits

    script:
    """
    echo -e "Fam\tHit\tProb\tE-value\tP-value\tScore\tSS\tCols\tQuery HMM\tTemplate HMM" > self_hits.tsv

    for hhr_file in ${hhr_folder}/*; do
        name=\$(basename \$hhr_file .hhr)
        awk '/^ No Hit/{flag=1;next}/^\$/{flag=0}flag' \$hhr_file | \
        awk -v name=\$name 'BEGIN {OFS="\t"} {val=substr(\$0, 42, 8) + 0; if (val < 0.001) { \$1=name; print \$0 }}' >> self_hits.tsv
    done
    """
}
