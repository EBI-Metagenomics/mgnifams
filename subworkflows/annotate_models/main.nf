#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT  } from "${params.moduleDir}/hhsuite/reformat/main.nf"
include { HHSUITE_HHBLITS   } from "${params.moduleDir}/hhsuite/hhblits/main.nf"
include { HHSUITE_HHSEARCH  } from "${params.moduleDir}/hhsuite/hhsearch/main.nf"
include { FILTER_HH_RESULTS } from "${params.moduleDir}/hhsuite/filter_hh_results.nf"

workflow ANNOTATE_MODELS {
    take:
    fasta_ch
    hh_mode
    
    main:
    a2m_fasta_ch = HHSUITE_REFORMAT(fasta_ch, "fas", "a3m").fa
    if (hh_mode == "hhblits") {
        hhr_ch = HHSUITE_HHBLITS(a2m_fasta_ch, params.hhdb_folder_path, params.db_name).hhr
    } else if (hh_mode == "hhsearch") {
        hhr_ch = HHSUITE_HHSEARCH(a2m_fasta_ch, params.hhdb_folder_path, params.db_name).hhr
    } else {
        throw new Exception("Invalid hh_mode value. Should be 'hhblits' or 'hhsearch'.")
    }
    FILTER_HH_RESULTS(hhr_ch)
}
