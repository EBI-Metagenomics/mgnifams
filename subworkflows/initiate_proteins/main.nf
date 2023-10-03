#!/usr/bin/env nextflow

include { KEEP_UNANNOTATED } from "$baseDir/modules/parsing.nf"
include { EXPORT_PROTEINS_CSV } from "$baseDir/modules/exporting.nf"

workflow initiate_proteins {
    main:
    def mgy90_path = params.mgy90_dataDir + params.mgy90_name

    Channel
        .fromPath(mgy90_path)
        .set { mgy90_file }

    fasta_file = KEEP_UNANNOTATED(mgy90_file)
    EXPORT_PROTEINS_CSV(fasta_file)

    emit:
    fasta_file
}