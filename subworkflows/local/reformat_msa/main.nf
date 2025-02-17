#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT    } from "../../../modules/local/hhsuite/reformat/main"
include { TRANSLATE_MSA_MGYPS } from "../../../modules/local/translate_msa_mgyps/main"

workflow REFORMAT_MSA {
    take:
    sto_ch
    
    main:
    ch_versions = Channel.empty()

    tmp_fa_ch = HHSUITE_REFORMAT(sto_ch, "sto", "fas").fa
    ch_versions = ch_versions.mix( HHSUITE_REFORMAT.out.versions )

    ch_fasta = TRANSLATE_MSA_MGYPS(tmp_fa_ch).fa
    // TODO ch_versions = ch_versions.mix( TRANSLATE_MSA_MGYPS.out.versions )

    emit:
    versions = ch_versions
    fasta    = ch_fasta
}
