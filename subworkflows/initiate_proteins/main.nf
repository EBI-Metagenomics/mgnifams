#!/usr/bin/env nextflow

include { FILTER_UNANNOTATED        } from "$launchDir/modules/initiate.nf"
include { FILTER_UNANNOTATED_SLICES } from "$launchDir/modules/initiate.nf"
include { EXPORT_PROTEINS_CSV       } from "$launchDir/modules/export.nf"

workflow INITIATE_PROTEINS {
    take:
    mode
    
    main:
    Channel
        .fromPath(params.mgy90_path)
        .set { mgy90_file }

    if (mode == "slice") {
        fasta = FILTER_UNANNOTATED_SLICES(mgy90_file, params.min_slice_length)
    } else {
        fasta = FILTER_UNANNOTATED(mgy90_file)
    }
    EXPORT_PROTEINS_CSV(fasta)

    emit:
    fasta
}
