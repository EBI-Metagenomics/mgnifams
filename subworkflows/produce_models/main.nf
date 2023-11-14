#!/usr/bin/env nextflow

include { MAFFT } from "$baseDir/modules/mafft.nf"
include { HMMBUILD } from "$baseDir/modules/hmm.nf"

workflow produce_models {
    take:
    familyFile

    main:
    mafft_folder = MAFFT(familyFile)
    build_folder = HMMBUILD(mafft_folder)

    emit:
    mafft_folder
    build_folder
}