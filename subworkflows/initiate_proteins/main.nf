#!/usr/bin/env nextflow

include { FILTER_UNANNOTATED    } from "$baseDir/modules/initiate.nf"
include { EXPORT_PROTEINS_CSV } from "$baseDir/modules/export.nf"

workflow INITIATE_PROTEINS {
    main:
    Channel
        .fromPath(mgy90_path)
        .set { mgy90_file }

    fasta_file = FILTER_UNANNOTATED(mgy90_file)
    EXPORT_PROTEINS_CSV(fasta_file)

    emit:
    fasta_file
}
