#!/usr/bin/env nextflow

include { PREPROCESS_INPUT   } from "./preprocess_input.nf"
include { INITIATE_PROTEINS  } from "./initiate_proteins.nf"
include { EXECUTE_CLUSTERING } from "./execute_clustering.nf"

workflow SETUP_CLUSTERS {
    take:
    input
    fasta_input_mode
    compress_mode

    main:
    if (!fasta_input_mode) {
        processed_input_protein_ch = PREPROCESS_INPUT(input, compress_mode).processed_input_protein_ch
        mgnifams_input_fa          = INITIATE_PROTEINS(processed_input_protein_ch).fasta_ch
    } else {
        mgnifams_input_fa = channel.fromPath(input)
    }
    clusters     = EXECUTE_CLUSTERING(mgnifams_input_fa)
    clusters_tsv = clusters.clusters_tsv.map { meta, filepath -> filepath }

    emit:
    mgnifams_input_fa
    clusters_tsv
}
