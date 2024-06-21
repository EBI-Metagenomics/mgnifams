#!/usr/bin/env nextflow

include { GENERATE_FAMILIES_PARALLEL } from "${projectDir}/../../subworkflows/generate_families_parallel/main.nf"
include { FLAG_TRANSMEMBRANE         } from "${projectDir}/subworkflows/flag_transmembrane/main.nf"
include { MOVE_TO_DIR                } from "${projectDir}/../../modules/family/main.nf"
include { REMOVE_REDUNDANCY          } from "${projectDir}/../../subworkflows/remove_redundancy/main.nf"

workflow {
    families_ch = GENERATE_FAMILIES_PARALLEL( params.clusters_tsv, params.checked_clusters_txt, params.mgnifams_input_fasta_path)

    msa_sto_ch = families_ch.msa_sto.collect()
    msa_sto_ch
        .map { files ->
            return [ [id:"flag_tm"], files ]
        }
        .set { tm_msa_ch }
    
    tm_ids_ch = FLAG_TRANSMEMBRANE(tm_msa_ch)
    
    seed_msa_sto_ch = families_ch.seed_msa_sto.collect()
    seed_msa_sto_dir = MOVE_TO_DIR(seed_msa_sto_ch, "seed_msa_sto")
    seed_msa_sto_dir
        .map { filepath ->
            return [ [id:"remove_redundancy"], file(filepath) ]
        }
        .set { seed_msa_sto_dir }
        
    generated_families = REMOVE_REDUNDANCY(seed_msa_sto_dir, seed_msa_sto_ch, \
        msa_sto_ch, families_ch.hmm.collect(), \
        families_ch.rf.collect(), families_ch.domtblout.collect(), families_ch.tsv.collect(), \
        families_ch.discarded.collect(), families_ch.successful.collect(), families_ch.converged.collect(), \
        families_ch.metadata.collect(), families_ch.logs.collect(), tm_ids_ch)
}
