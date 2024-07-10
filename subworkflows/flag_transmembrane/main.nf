#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES } from "../../modules/family/main.nf"
include { DEEPTMHMM                         } from "../../modules/deeptmhmm/main.nf"
include { FLAG_TM                           } from "../../modules/deeptmhmm/flag_tm.nf"

workflow FLAG_TRANSMEMBRANE {
    take:
    tm_msa_ch

    main:
    fa_ch     = EXTRACT_FIRST_STOCKHOLM_SEQUENCES(tm_msa_ch)
    gff3_ch   = DEEPTMHMM(fa_ch).gff3
    tm_ids_ch = FLAG_TM(gff3_ch, params.tm_fraction_threshold)

    emit:
    fa_ch     = fa_ch
    tm_ids_ch = tm_ids_ch
}
