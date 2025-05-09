#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT  } from "../modules/nf-core/hhsuite/reformat/main.nf"
include { HHSUITE_HHBLITS   } from "../modules/nf-core/hhsuite/hhblits/main.nf"
include { HHSUITE_HHSEARCH  } from "../modules/nf-core/hhsuite/hhsearch/main.nf"
include { FILTER_HH_RESULTS } from "../modules/nf-core/hhsuite/filter_hh_results.nf"

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
