#!/usr/bin/env nextflow

include { PIGZ_UNCOMPRESS            } from "../../../modules/nf-core/pigz/uncompress/main"
include { BZ2_UNCOMPRESS             } from "../../../modules/local/bz2_uncompress/main"
include { EXTRACT_UNANNOTATED_SLICES } from "../../../modules/local/extract_unannotated_slices/main"

workflow EXTRACT_UNANNOTATED_FASTA {
    take:
    sequence_explorer_protein_ch
    compress_mode
    
    main:
    if (compress_mode == 'gz') {
        sequence_explorer_protein_ch = PIGZ_UNCOMPRESS( sequence_explorer_protein_ch ).file
    } else if (compress_mode == 'bz2') {
        sequence_explorer_protein_ch = BZ2_UNCOMPRESS( sequence_explorer_protein_ch ).file
    }

    sequence_explorer_protein_ch
        .splitText(file:true, by: params.input_csv_chunk_size, keepHeader: true)
        .map { meta, file ->
            [[id: meta.id, chunk: file.getBaseName(1).split('\\.')[-1]], file]
        }
        .set { sequence_chunk_ch }
    
    fasta_chunk_ch = EXTRACT_UNANNOTATED_SLICES( sequence_chunk_ch, params.min_sequence_length ).fa
    
    fasta_chunk_ch
        .map { meta, file ->
            file
        }
        .collectFile(name: "mgnifams_input.fa", storeDir: params.outdir)
        .map { file ->
            [[id: 'mgnifams_input'], file]
        }
        .set { fasta_ch }

    emit:
    fasta_ch
}
