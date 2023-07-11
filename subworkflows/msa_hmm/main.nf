#!/usr/bin/env nextflow

include { MAFFT } from "$baseDir/modules/mafft/mafft.nf"
include { HMMBUILD; HMMSCAN } from "$baseDir/modules/hmm/hmm.nf"
include { SLICE } from "$baseDir/modules/hmm/slice.nf"

workflow msa_hmm {
    take:
    familyFile
    uniprot_sprot_fasta_path

    main:
    mafft_ch = MAFFT(familyFile)
    build_ch = HMMBUILD(mafft_ch)
    tblout_ch = HMMSCAN(build_ch, uniprot_sprot_fasta_path.first()).tblout_ch
    SLICE(tblout_ch)
}