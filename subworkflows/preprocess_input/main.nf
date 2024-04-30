#!/usr/bin/env nextflow

include { DECOMPRESS_GZ  } from "${params.moduleDir}/preprocess.nf"
include { DECOMPRESS_BZ2 } from "${params.moduleDir}/preprocess.nf"
include { REMOVE_HEADER  } from "${params.moduleDir}/preprocess.nf"

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
    preprocessed_sequence_explorer_protein_ch = REMOVE_HEADER(sequence_explorer_protein_ch)

    emit:
    preprocessed_sequence_explorer_protein_ch
}
