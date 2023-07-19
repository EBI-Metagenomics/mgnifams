#!/usr/bin/env nextflow

include { MAFFT } from "$baseDir/modules/mafft.nf"
include { HMMBUILD } from "$baseDir/modules/hmm.nf"

workflow produce_models {
    take:
    familyFile

    main:
    mafft_ch = MAFFT(familyFile)
    build_ch = HMMBUILD(mafft_ch)

    emit:
    mafft_ch
    build_ch
}