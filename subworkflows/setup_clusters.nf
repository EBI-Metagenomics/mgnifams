#!/usr/bin/env nextflow

include { EXTRACT_UNANNOTATED_FASTA } from "../subworkflows/extract_unannotated_fasta.nf"
include { EXECUTE_CLUSTERING        } from "../subworkflows/execute_clustering.nf"

workflow SETUP_CLUSTERS {
    take:
    input
    fasta_input_mode
    compress_mode

    main:
    if (!fasta_input_mode) {
        mgnifams_input_fa = EXTRACT_UNANNOTATED_FASTA( input, compress_mode ).fasta_ch
    } else {
        mgnifams_input_fa = channel.fromPath(input)
    }

    // TODO
    mgnifams_input_fa.view()
    // clusters     = EXECUTE_CLUSTERING(mgnifams_input_fa)
    // clusters_tsv = clusters.clusters_tsv.map { meta, filepath -> filepath }

    // emit:
    // mgnifams_input_fa
    // clusters_tsv
}
