include { PIGZ_UNCOMPRESS            } from "../../../modules/nf-core/pigz/uncompress/main"
include { BZ2_UNCOMPRESS             } from "../../../modules/local/bz2_uncompress/main"
include { EXTRACT_UNANNOTATED_SLICES } from "../../../modules/local/extract_unannotated_slices/main"

workflow EXTRACT_UNANNOTATED_FASTA {
    take:
    ch_sequence_explorer_protein
    compress_mode
    input_csv_chunk_size
    min_sequence_length
    outdir
    
    main:
    ch_versions = Channel.empty()

    if (compress_mode == 'gz') {
        ch_sequence_explorer_protein = PIGZ_UNCOMPRESS( ch_sequence_explorer_protein ).file
        ch_versions = ch_versions.mix( PIGZ_UNCOMPRESS.out.versions )
    } else if (compress_mode == 'bz2') {
        ch_sequence_explorer_protein = BZ2_UNCOMPRESS( ch_sequence_explorer_protein ).file
        ch_versions = ch_versions.mix( BZ2_UNCOMPRESS.out.versions )
    }

    ch_sequence_explorer_chunk = ch_sequence_explorer_protein
        .splitText(file:true, by: input_csv_chunk_size, keepHeader: true)
        .map { meta, file ->
            [[id: meta.id, chunk: file.getBaseName(1).split('\\.')[-1]], file]
        }
    
    ch_fasta_chunk = EXTRACT_UNANNOTATED_SLICES( ch_sequence_explorer_chunk, min_sequence_length ).fa
    ch_versions = ch_versions.mix( EXTRACT_UNANNOTATED_SLICES.out.versions )
    
    ch_fasta = ch_fasta_chunk
        .map { meta, file ->
            file
        }
        .collectFile(name: "mgnifams_v2.fa", storeDir: outdir)
        .map { file ->
            [[id: 'mgnifams_v2'], file]
        }

    emit:
    versions = ch_versions
    fasta    = ch_fasta
}
