#!/bin/bash

# filter_clusters.sh
# Usage: filter_clusters.sh [input_tsv] [threshold] [output_tsv]

input_tsv="$1"
threshold="$2"
output_tsv="$3"

awk -F'\t' '{
    if ($1 != prev && NR > 1) {
        if (count >= threshold) 
            printf "%s", lines;
        lines = "";
        count = 0;
    }
    lines = lines $0 "\n";
    count++;
    prev = $1;
} 
END {
    if (count >= threshold) 
        printf "%s", lines;
}' threshold="$threshold" "$input_tsv" > "$output_tsv"
