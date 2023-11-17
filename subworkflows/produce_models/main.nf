#!/usr/bin/env nextflow

include { MAFFT          } from "../../modules/mafft/main.nf"
include { HMMER_HMMBUILD } from "../../modules/hmmer/hmmbuild/main.nf"

workflow PRODUCE_MODELS {
    take:
    familyFile

    main:
    input = familyFile.map { file ->
        return [ [id:file.getBaseName()], file]
    }
    fas_ch = MAFFT( input, [] ).fas
    hmm_ch = HMMER_HMMBUILD( fas_ch, [] ).hmm

    emit:
    fas_ch
    hmm_ch
}