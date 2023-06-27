#!/usr/bin/env nextflow

include { MAFFT } from '../../modules/mafft/mafft.nf'
include { HMMBUILD } from '../../modules/hmm/hmm.nf'

workflow msa_hmm {
    take: familyFile
    main:
        MAFFT(familyFile) | HMMBUILD 
}