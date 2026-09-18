/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UPDATE_MGNIFAMS: refresh full MSAs and representatives of existing families
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UPDATE_MGNIFAMS {

    take:
    ch_samplesheet // channel: [ meta, sequences_parquet, pfam_parquet, hmm_lib, db_config or [] ]

    main:
    ch_samplesheet.view { meta, sequences, pfam, hmms, _db_config -> "update_mgnifams ${meta.id}: ${sequences.name}, ${pfam.name}, ${hmms.name}" }

    emit:
    channel.empty() // multiqc_report
}
