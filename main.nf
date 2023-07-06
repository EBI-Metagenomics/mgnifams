#!/usr/bin/env nextflow

include { mmseqs2 } from "$baseDir/subworkflows/mmseqs2/main.nf"
include { fasta_families } from "$baseDir/subworkflows/fasta_families/main.nf"
include { msa_hmm } from "$baseDir/subworkflows/msa_hmm/main.nf"

workflow {
    Channel
        .fromPath(params.fasta_path) 
        .set { fastaFile }

    Channel
        .fromPath(params.uniprot_sprot_fasta_path) 
        .set { uniprot_sprot_fasta_path }

    cluster_tsv_ch = mmseqs2(fastaFile)
    families_ch = fasta_families(cluster_tsv_ch, fastaFile)
    msa_hmm(families_ch, uniprot_sprot_fasta_path)
}