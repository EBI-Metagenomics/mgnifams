#!/usr/bin/env nextflow

include { HMMSCAN } from "$baseDir/modules/hmm.nf"
include { INTERPROSCAN } from "$baseDir/modules/interproscan.nf"
include { EGGNOG_MAPPER } from "$baseDir/modules/eggnog.nf"
// TODO protENN

process getTopDebugLines { // TODO remove in Codon
    input:
    path(fasta)

    output:
    path("top_reps.fasta")

    script:
    lines = params.debug_top_n_reps * 2
    """
    head -n ${lines} ${fasta} > top_reps.fasta
    """
}

workflow annotate_unknown {
    take:
    unknown_reps_fasta
    build_ch
    
    main:
    Channel
        .fromPath(params.uniprot_sprot_fasta_path) 
        .set { uniprot_sprot_fasta_path }

    top_reps = getTopDebugLines(unknown_reps_fasta) // TODO remove
    INTERPROSCAN(top_reps) // TODO: reps_fa
    EGGNOG_MAPPER(top_reps) // TODO:  reps_fa
    tblout_ch = HMMSCAN(build_ch, uniprot_sprot_fasta_path.first()).tblout_ch

    emit:
    tblout_ch
}