#!/usr/bin/env nextflow

include { SLICE; PRINT_SLICED; CONCAT_FASTA } from "$baseDir/modules/slicing/slicing.nf"

workflow slice_unannotated {
    take:
    msa_ch
    tblout_ch
    
    main:
    // Map the msa_ch channel to include the key (base name)
    msa_mapped_ch = msa_ch.map { file -> [ file.baseName.split('_')[0], file ] }

    slice_ch = SLICE(tblout_ch)

    // Map the slice_ch channel to include the key (base name)
    slice_mapped_ch = slice_ch.map { file -> [ file.baseName.split('_')[0], file ] }

    // Join the mapped channels based on the key
    matched_files_ch = msa_mapped_ch.join(slice_mapped_ch)

    // Use the joined channel in your processes, counts as 1 input, tuple
    sliced_ch = PRINT_SLICED(matched_files_ch.map { it -> [ it[1], it[2] ] }) // Pass msa_file and sliced_file to the PRINT_SLICED process, it[0] is the key
    fasta_file = CONCAT_FASTA(sliced_ch.collect())

    emit:
    fasta_file
}