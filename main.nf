#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ebi-metagenomics/mgnifams
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ebi-metagenomics/mgnifams
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MGNIFAMS                } from './workflows/mgnifams'
include { INIT_DB                 } from './workflows/init_db'
include { UPDATE_DB               } from './workflows/update_db'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_mgnifams_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_mgnifams_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow EBIMETAGENOMICS_MGNIFAMS {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_multiqc = Channel.empty()

    //
    // WORKFLOW: Run main pipeline
    //
    if (params.mode == 'run_mgnifams_pipeline') {
        MGNIFAMS( 
            samplesheet, params.fasta_input_mode, params.compress_mode, \
            params.input_csv_chunk_size, params.min_sequence_length, params.outdir, \
            params.minimum_members, params.clusters_chunk_size, params.mgnifams_discard_min_rep_length, \
            params.mgnifams_discard_max_rep_length, params.mgnifams_discard_min_starting_membership, \
            params.mgnifams_max_seq_identity, params.mgnifams_max_seed_seqs, params.mgnifams_max_gap_occupancy, \
            params.mgnifams_recruit_evalue_cutoff, params.mgnifams_recruit_hit_length_percentage, \
            params.redundant_length_threshold, params.redundant_score_threshold, \
            params.similarity_score_threshold, params.starting_id, \
            params.pdb_chunk_size, params.esmfold_db, params.esmfold_params_path, \
            params.esmfold_3B_v1, params.esm2_t36_3B_UR50D, params.esm2_t36_3B_UR50D_contact_regression, \
            params.num_recycles_esmfold, params.pdb_chunk_size_long, \
            params.skip_deeptmhmm, params.deeptmhmm_path, params.pfam_path, params.funfams_path, \
            params.hh_mode, params.hhdb_path, params.foldseek_db_path, params.query_hmm_length_threshold, \
            params.multiqc_config, params.multiqc_logo, params.multiqc_methods_description
        )
        ch_multiqc = MGNIFAMS.out.multiqc_report
    }
    //
    // WORKFLOW: Run initialize mgnifams db workflow
    //
    else if (params.mode == 'init_mgnifams_db') {
        INIT_DB( 
            samplesheet
        )
    }
    //
    // WORKFLOW: Run update mgnifams db workflow
    //
    else if (params.mode == 'update_mgnifams_db') {
        UPDATE_DB( 
            samplesheet, params.query_result_chunks
        )
    }

    emit:
    multiqc_report = ch_multiqc // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //

    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )
    
    //
    // WORKFLOW: Run main workflow
    //
    EBIMETAGENOMICS_MGNIFAMS (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        EBIMETAGENOMICS_MGNIFAMS.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/