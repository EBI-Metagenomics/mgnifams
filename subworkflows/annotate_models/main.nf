#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT } from "${launchDir}/modules/hhsuite/reformat/main.nf"
include { HHSUITE_HHBLITS  } from "${launchDir}/modules/hhsuite/hhblits/main.nf"

workflow ANNOTATE_MODELS {
    take:
    fasta_ch
    db
    
    main:
    a2m_fasta_ch = HHSUITE_REFORMAT(fasta_ch, "sto", "a3m").fa
    hhr_ch = HHSUITE_HHBLITS(a2m_fasta_ch, db, "pfam").hhr // TODO check if .first is needed in db

    emit:
    hhr_ch
}