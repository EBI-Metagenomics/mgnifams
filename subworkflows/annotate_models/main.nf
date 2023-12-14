#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT  } from "${launchDir}/modules/hhsuite/reformat/main.nf"
include { HHSUITE_HHBLITS   } from "${launchDir}/modules/hhsuite/hhblits/main.nf"
include { HHSUITE_HHSEARCH  } from "${launchDir}/modules/hhsuite/hhsearch/main.nf"
include { FILTER_HH_RESULTS } from "${launchDir}/modules/hhsuite/filter.nf"

workflow ANNOTATE_MODELS {
    take:
    fasta_ch
    db
    mode
    
    main:
    a2m_fasta_ch = HHSUITE_REFORMAT(fasta_ch, "fas", "a3m").fa
    if (mode == "hhblits") {
        hhr_ch = HHSUITE_HHBLITS(a2m_fasta_ch, db, "pfam").hhr
    } else if (mode == "hhsearch") {
        hhr_ch = HHSUITE_HHSEARCH(a2m_fasta_ch, db, "pfam").hhr
    } else {
        throw new Exception("Invalid mode value. Should be 'hhblits' or 'hhsearch'.")
    }
    results = FILTER_HH_RESULTS(hhr_ch)

    emit:
    unannotated = results.unannotated
}
