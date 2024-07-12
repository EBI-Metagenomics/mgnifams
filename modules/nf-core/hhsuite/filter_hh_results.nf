process FILTER_HH_RESULTS {
    label "general"

    input:
    tuple val(meta), path(hhr_folder)

    output:
    tuple val(meta), path("pfam_hits"), optional: true, emit: pfam_hits
    path("annotated_models.txt")      , optional: true, emit: annotated_models

    script:
    """
    mkdir -p pfam_hits

    for hhr_file in ${hhr_folder}/*; do
        name=\$(basename \$hhr_file .hhr)
        awk '/^ No Hit/{flag=1;next}/^\$/{flag=0}flag' \$hhr_file | awk '{val=substr(\$0, 42, 8) + 0; if (val < 0.001) print \$0}' > pfam_hits/\$name.txt
        if [ -s pfam_hits/\$name.txt ]; then
            echo \$name >> annotated_models.txt
        else
            rm pfam_hits/\$name.txt
        fi
    done
    """
}
