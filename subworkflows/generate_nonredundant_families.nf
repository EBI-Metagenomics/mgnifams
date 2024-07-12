#!/usr/bin/env nextflow

include { GENERATE_FAMILIES_PARALLEL } from "./generate_families_parallel.nf"
include { FLAG_TRANSMEMBRANE         } from "./flag_transmembrane.nf"
include { MOVE_TO_DIR                } from "../modules/local/move_to_dir.nf"
include { REMOVE_REDUNDANCY          } from "./remove_redundancy.nf"

workflow GENERATE_NONREDUNDANT_FAMILIES {
    take:
    clusters_tsv
    checked_clusters
    mgnifams_input_fa

    main:
    families_ch = GENERATE_FAMILIES_PARALLEL(clusters_tsv, [], mgnifams_input_fa)

    msa_sto_ch = families_ch.msa_sto.collect()
    msa_sto_ch
        .map { files ->
            return [ [id:"flag_tm"], files ]
        }
        .set { tm_msa_ch }
    
    tm_ch     = FLAG_TRANSMEMBRANE(tm_msa_ch)
    rep_fa_ch = tm_ch.fa_ch
    tm_ids_ch = tm_ch.tm_ids_ch
    
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
        families_ch.metadata.collect(), families_ch.logs.collect(), tm_ids_ch, rep_fa_ch)

    emit:
    seed_msa_sto = generated_families.seed_msa_sto
    msa_sto      = generated_families.msa_sto
    metadata     = generated_families.metadata
    converged    = generated_families.converged
    tsv          = generated_families.tsv
    hmm          = generated_families.hmm
    rf           = generated_families.rf
}
