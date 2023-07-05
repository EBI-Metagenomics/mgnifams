#!/usr/bin/env nextflow

params.fasta_path = "$baseDir/data/input/mgnify_500K_proteins.fa.gz"

include { mmseqs2 } from "$baseDir/subworkflows/mmseqs2/main.nf"
include { fasta_families } from "$baseDir/subworkflows/fasta_families/main.nf"
include { msa_hmm } from "$baseDir/subworkflows/msa_hmm/main.nf"

workflow {
    Channel
        .fromPath(params.fasta_path) 
        .set { fastaFile }

    cluster_tsv_ch = mmseqs2(fastaFile)
    families_ch = fasta_families(cluster_tsv_ch, fastaFile)
    msa_hmm(families_ch)
}