#!/usr/bin/env nextflow

include { FIND_FASTA_BY_ID } from "$baseDir/modules/general.nf"
include { ESMFOLD } from "$baseDir/modules/esmfold/main.nf"

workflow annotate_structures {
    take:
    ids
    fasta
    
    main:
    FIND_FASTA_BY_ID(ids, fasta)

    emit:
    FIND_FASTA_BY_ID.out
}