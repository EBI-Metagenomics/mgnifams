#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT    } from "${params.moduleDir}/hhsuite/reformat/main.nf"
include { TRANSLATE_MSA_MGYPS } from "${params.moduleDir}/translate.nf"

workflow REFORMAT_MSA {
    take:
    sto_ch
    
    main:
    tmp_fa_ch = HHSUITE_REFORMAT(sto_ch, "sto", "fas").fa
    fa_ch = TRANSLATE_MSA_MGYPS(tmp_fa_ch).fa

    emit:
    fa_ch
}
