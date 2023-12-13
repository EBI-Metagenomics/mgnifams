process FILTER_HH_RESULTS {
    publishDir "${params.outdir}/hh", mode: "copy"
    label "general"

    input:
    tuple val(meta), path(hhr_folder)

    output:
    tuple val(meta), path("hits"), optional: true, emit: hits
    path("unannotated.txt")      , optional: true, emit: unannotated

    script:
    """
    mkdir -p hits

    for hhr_file in ${hhr_folder}/*; do
        name=\$(basename \$hhr_file .hhr)
        awk '/^ No Hit/{flag=1;next}/^\$/{flag=0}flag' \$hhr_file | awk '{val=substr(\$0, 42, 8) + 0; if (val < 0.001) print \$0}' > hits/\$name.tsv
        if [ ! -s hits/\$name.tsv ]; then
            echo \$name >> unannotated.txt
            rm hits/\$name.tsv
        fi
    done
    """
}
