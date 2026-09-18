/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UPDATE_MGNIFAMS: refresh full MSAs and representatives of existing families
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UPDATE_FAMILIES } from '../subworkflows/local/update_families'

workflow UPDATE_MGNIFAMS {

    take:
    ch_samplesheet      // channel: [ meta, sequences_parquet, pfam_parquet, hmm_lib, db_config or [] ]
    parquet_chunks
    min_sequence_length
    hmm_chunk_size
    outdir

    main:
    UPDATE_FAMILIES( ch_samplesheet, parquet_chunks, min_sequence_length, hmm_chunk_size, outdir )

    emit:
    channel.empty() // multiqc_report
}
