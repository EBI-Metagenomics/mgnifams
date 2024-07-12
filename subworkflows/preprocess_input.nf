#!/usr/bin/env nextflow

include { DECOMPRESS_GZ  } from "../modules/local/decompress_gz.nf"
include { DECOMPRESS_BZ2 } from "../modules/local/decompress_bz2.nf"
include { REMOVE_HEADER  } from "../modules/local/remove_header.nf"

workflow PREPROCESS_INPUT {
    take:
    sequence_explorer_protein_ch
    compress_mode
    
    main:
    if (compress_mode == 'gz') {
        sequence_explorer_protein_ch = DECOMPRESS_GZ(sequence_explorer_protein_ch)
    } else if (compress_mode == 'bz2') {
        sequence_explorer_protein_ch = DECOMPRESS_BZ2(sequence_explorer_protein_ch)
    }
    processed_input_protein_ch = REMOVE_HEADER(sequence_explorer_protein_ch)

    emit:
    processed_input_protein_ch
}
