#!/usr/bin/env nextflow

include { CONCAT_FILES } from "$baseDir/modules/general.nf"
include { EXPORT_PROTEINS_CSV } from "$baseDir/modules/exporting.nf"

workflow initiate_proteins {
    main:
    def mgnify_fasta_path = params.dataDir + params.mgnify_fasta_name
    def uniprot_sprot_fasta_input_debug_path = params.dataDir + params.uniprot_sprot_fasta_input_debug_name

    Channel
        .fromPath(mgnify_fasta_path) 
        .set { mgnify_fasta_file }
    Channel
        .fromPath(uniprot_sprot_fasta_input_debug_path) 
        .set { uniprot_sprot_fasta_input_debug_file } // TODO change to main Uniprot file

    combined_fasta_file = CONCAT_FILES(mgnify_fasta_file, uniprot_sprot_fasta_input_debug_file)
    EXPORT_PROTEINS_CSV(combined_fasta_file)

    emit:
    combined_fasta_file
}