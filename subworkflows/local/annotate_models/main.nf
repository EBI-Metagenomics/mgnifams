#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT  } from "../../../modules/nf-core/hhsuite/reformat/main"
include { HHSUITE_HHBLITS   } from "../../../modules/local/hhsuite/hhblits/main"
include { HHSUITE_HHSEARCH  } from "../../../modules/local/hhsuite/hhsearch/main"
include { FILTER_HH_RESULTS } from "../../../modules/local/hhsuite/filter_hh_results/main"

workflow ANNOTATE_MODELS {
    take:
    seed_msa
    hh_mode
    hhdb_path
    hhdb_name
    
    main:
    ch_versions = Channel.empty()

    HHSUITE_REFORMAT( seed_msa, "fas", "a3m" )
    ch_versions = ch_versions.mix( HHSUITE_REFORMAT.out.versions )

    if (hh_mode == "hhblits") {
        hhr_ch = HHSUITE_HHBLITS( HHSUITE_REFORMAT.out.msa, hhdb_path, hhdb_name ).hhr
        ch_versions = ch_versions.mix( HHSUITE_HHBLITS.out.versions )
    } else if (hh_mode == "hhsearch") {
        hhr_ch = HHSUITE_HHSEARCH( HHSUITE_REFORMAT.out.msa, hhdb_path, hhdb_name ).hhr
        ch_versions = ch_versions.mix( HHSUITE_HHSEARCH.out.versions )
    } else {
        throw new Exception("Invalid hh_mode value. Should be 'hhblits' or 'hhsearch'.")
    }
    // ch_pfam_hits = FILTER_HH_RESULTS(hhr_ch).pfam_hits
    // TODO ch_versions = ch_versions.mix( FILTER_HH_RESULTS.out.versions )

    emit:
    versions  = ch_versions
    // pfam_hits = ch_pfam_hits
}
