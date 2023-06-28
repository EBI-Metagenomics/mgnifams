#!/usr/bin/env nextflow

include { MAFFT } from "$baseDir/modules/mafft/mafft.nf"
include { HMMBUILD } from "$baseDir/modules/hmm/hmm.nf"

workflow msa_hmm {
    take: familyFile
    main:
        MAFFT(familyFile) | HMMBUILD 
}