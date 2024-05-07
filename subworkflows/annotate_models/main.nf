#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT  } from "${params.moduleDir}/hhsuite/reformat/main.nf"
include { HHSUITE_HHBLITS   } from "${params.moduleDir}/hhsuite/hhblits/main.nf"
include { HHSUITE_HHSEARCH  } from "${params.moduleDir}/hhsuite/hhsearch/main.nf"
include { FILTER_HH_RESULTS } from "${params.moduleDir}/hhsuite/filter_hh_results.nf"

workflow ANNOTATE_MODELS {
    take:
    fasta_ch
    mode
    
    main:
    a2m_fasta_ch = HHSUITE_REFORMAT(fasta_ch, "fas", "a3m").fa
    if (mode == "hhblits") {
        hhr_ch = HHSUITE_HHBLITS(a2m_fasta_ch, params.hhdb_folder_path, "pfam").hhr
    } else if (mode == "hhsearch") {
        hhr_ch = HHSUITE_HHSEARCH(a2m_fasta_ch, params.hhdb_folder_path, "pfam").hhr
    } else {
        throw new Exception("Invalid mode value. Should be 'hhblits' or 'hhsearch'.")
    }
    results = FILTER_HH_RESULTS(hhr_ch)

    emit:
    unannotated = results.unannotated
}
