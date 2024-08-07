#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES } from "${params.moduleDir}/family/main.nf"
include { DEEPTMHMM                         } from "${params.moduleDir}/deeptmhmm/main.nf"
include { FLAG_TM                           } from "${params.moduleDir}/deeptmhmm/flag_tm.nf"

workflow FLAG_TRANSMEMBRANE {
    take:
    tm_msa_ch

    main:
    extracted_res = EXTRACT_FIRST_STOCKHOLM_SEQUENCES(tm_msa_ch)
    fa_ch         = extracted_res.fa
    prob_ids      = extracted_res.prob_ids
    gff3_ch       = DEEPTMHMM(fa_ch).gff3
    tm_ids_ch     = FLAG_TM(gff3_ch, params.tm_fraction_threshold)

    emit:
    fa_ch     = fa_ch
    tm_ids_ch = tm_ids_ch
    prob_ids  = prob_ids
}
