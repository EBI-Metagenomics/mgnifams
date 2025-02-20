#!/usr/bin/env nextflow

include { PIGZ_UNCOMPRESS            } from "../../../modules/nf-core/pigz/uncompress/main"
include { BZ2_UNCOMPRESS             } from "../../../modules/local/bz2_uncompress/main"
include { EXTRACT_UNANNOTATED_SLICES } from "../../../modules/local/extract_unannotated_slices/main"

workflow EXTRACT_UNANNOTATED_FASTA {
    take:
    ch_sequence_explorer_protein
    
    main:
    ch_versions = Channel.empty()

    if (params.compress_mode == 'gz') {
        ch_sequence_explorer_protein = PIGZ_UNCOMPRESS( ch_sequence_explorer_protein ).file
        ch_versions = ch_versions.mix( PIGZ_UNCOMPRESS.out.versions )
    } else if (params.compress_mode == 'bz2') {
        ch_sequence_explorer_protein = BZ2_UNCOMPRESS( ch_sequence_explorer_protein ).file
        ch_versions = ch_versions.mix( BZ2_UNCOMPRESS.out.versions )
    }

    ch_sequence_explorer_chunk = ch_sequence_explorer_protein
        .splitText(file:true, by: params.input_csv_chunk_size, keepHeader: true)
        .map { meta, file ->
            [[id: meta.id, chunk: file.getBaseName(1).split('\\.')[-1]], file]
        }
    
    ch_fasta_chunk = EXTRACT_UNANNOTATED_SLICES( ch_sequence_explorer_chunk, params.min_sequence_length ).fa
    ch_versions = ch_versions.mix( EXTRACT_UNANNOTATED_SLICES.out.versions )
    
    ch_fasta = ch_fasta_chunk
        .map { meta, file ->
            file
        }
        .collectFile(name: "mgnifams_input.fa", storeDir: params.outdir)
        .map { file ->
            [[id: 'mgnifams_input'], file]
        }

    emit:
    versions = ch_versions
    fasta    = ch_fasta
}
