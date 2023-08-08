#!/usr/bin/env nextflow

include { CONCAT_FILES } from "$baseDir/modules/general.nf"
include { EXPORT_PROTEINS_CSV } from "$baseDir/modules/exporting.nf"

workflow initiate_proteins {
    main:
    def mgy_fasta_path = params.mgy_dataDir + params.mgy_fasta_name
    def uniprot_sprot_fasta_input_path = params.dataDir + params.uniprot_sprot_fasta_input_name

    Channel
        .fromPath(mgy_fasta_path)
        .collect()
        .set { mgy_folder }
    Channel
        .fromPath(uniprot_sprot_fasta_input_path) 
        .set { uniprot_sp }

    combined_fasta_file = CONCAT_FILES(mgy_folder, uniprot_sp)
    EXPORT_PROTEINS_CSV(combined_fasta_file)

    emit:
    combined_fasta_file
}