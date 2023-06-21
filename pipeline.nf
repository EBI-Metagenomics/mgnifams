#!/usr/bin/env nextflow

params.fasta_path = "data/mgnify_500K_proteins.fa.gz" // "data/cluster_reps_50k.fa"
params.family_path = "results/family.fasta"

include { CREATE_DB; LINCLUST; CREATE_TSV } from './modules/mmseqs2/mmseqs2.nf'
include { MAFFT } from './modules/mafft/mafft.nf'
include { HMMBUILD } from './modules/hmm/hmm.nf'

// workflow {
//     Channel
//         .fromPath(params.fasta_path) 
//         .set { fastaFile }

//     db_ch = CREATE_DB(fastaFile)
//     clust_ch = LINCLUST(db_ch)
//     CREATE_TSV(db_ch, clust_ch)
// }

workflow {
    Channel
        .fromPath(params.family_path) 
        .set { familyFile }

        MAFFT(familyFile) | HMMBUILD
        
}