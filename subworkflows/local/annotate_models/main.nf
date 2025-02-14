#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT  } from "../../../modules/local/hhsuite/reformat/main"
include { HHSUITE_HHBLITS   } from "../../../modules/local/hhsuite/hhblits/main"
include { HHSUITE_HHSEARCH  } from "../../../modules/local/hhsuite/hhsearch/main"
include { FILTER_HH_RESULTS } from "../../../modules/local/hhsuite/filter_hh_results/main"

workflow ANNOTATE_MODELS {
    take:
    fasta_ch
    
    main:
    a2m_fasta_ch = HHSUITE_REFORMAT(fasta_ch, "fas", "a3m").fa
    if (params.hh_mode == "hhblits") {
        hhr_ch = HHSUITE_HHBLITS(a2m_fasta_ch, params.hhdb_folder_path, params.db_name).hhr
    } else if (params.hh_mode == "hhsearch") {
        hhr_ch = HHSUITE_HHSEARCH(a2m_fasta_ch, params.hhdb_folder_path, params.db_name).hhr
    } else {
        throw new Exception("Invalid hh_mode value. Should be 'hhblits' or 'hhsearch'.")
    }
    pfam_hits = FILTER_HH_RESULTS(hhr_ch).pfam_hits

    emit:
    pfam_hits
}
