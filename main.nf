#!/usr/bin/env nextflow

params.fasta_path = "$baseDir/data/input/mgnify_500K_proteins.fa.gz" // "data/cluster_reps_50k.fa"
params.family_path = "$baseDir/data/output/family.fasta"

include { mmseqs2 } from "$baseDir/subworkflows/mmseqs2/main.nf"
include { msa_hmm } from "$baseDir/subworkflows/msa_hmm/main.nf"

workflow {
    Channel
        .fromPath(params.fasta_path) 
        .set { fastaFile }
    Channel
        .fromPath(params.family_path) 
        .set { familyFile }

    mmseqs2(fastaFile)
    msa_hmm(familyFile)
}