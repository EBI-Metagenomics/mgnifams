#!/usr/bin/env nextflow

include { PREPROCESS_INPUT   } from "./preprocess_input/main.nf"
include { INITIATE_PROTEINS  } from "./initiate_proteins/main.nf"
include { EXECUTE_CLUSTERING } from "./execute_clustering/main.nf"

workflow SETUP_CLUSTERS {
    take:
    sequence_explorer_protein_ch
    compress_mode

    main:
    processed_input_protein_ch = PREPROCESS_INPUT(sequence_explorer_protein_ch, compress_mode).processed_input_protein_ch
    mgnifams_input_fa          = INITIATE_PROTEINS(processed_input_protein_ch).fasta_ch
    clusters                   = EXECUTE_CLUSTERING(mgnifams_input_fa)
    clusters_tsv               = clusters.clusters_tsv.map { meta, filepath -> filepath }

    emit:
    mgnifams_input_fa
    clusters_tsv
}
