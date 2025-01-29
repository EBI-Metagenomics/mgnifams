#!/usr/bin/env nextflow

include { GENERATE_FAMILIES_PARALLEL } from "../subworkflows/generate_families_parallel.nf"
include { FLAG_TRANSMEMBRANE         } from "../subworkflows/flag_transmembrane.nf"
include { MOVE_TO_DIR                } from "../modules/local/move_to_dir.nf"
include { REMOVE_REDUNDANCY          } from "../subworkflows/remove_redundancy.nf"

workflow GENERATE_NONREDUNDANT_FAMILIES {
    take:
    clusters_tsv
    checked_clusters
    mgnifams_input_fa

    main:
    families_ch = GENERATE_FAMILIES_PARALLEL(clusters_tsv, [], mgnifams_input_fa)

    families_ch.msa_sto
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"msa_sto"], file ]
        }
        .set { msa_sto_ch }
    
    tm_ch     = FLAG_TRANSMEMBRANE(msa_sto_ch)
    rep_fa_ch = tm_ch.fa_ch
    tm_ids_ch = tm_ch.tm_ids_ch
    prob_ids  = tm_ch.prob_ids
    
    families_ch.seed_msa_sto
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"seed_msa_sto"], file ]
        }
        .set { seed_msa_sto_ch }

    seed_msa_sto_dir = MOVE_TO_DIR(seed_msa_sto_ch, "seed_msa_sto")

    // TODO change collects to .map because they name have meta
    // generated_families = REMOVE_REDUNDANCY(seed_msa_sto_dir, seed_msa_sto_ch, \
    //     msa_sto_ch, families_ch.hmm.collect(), \
    //     families_ch.rf.collect(), families_ch.domtblout.collect(), families_ch.tsv.collect(), \
    //     families_ch.discarded.collect(), families_ch.successful.collect(), families_ch.converged.collect(), \
    //     families_ch.metadata.collect(), families_ch.logs.collect(), tm_ids_ch, prob_ids, rep_fa_ch)

    // emit:
    // seed_msa_sto = generated_families.seed_msa_sto
    // msa_sto      = generated_families.msa_sto
    // metadata     = generated_families.metadata
    // converged    = generated_families.converged
    // tsv          = generated_families.tsv
    // hmm          = generated_families.hmm
    // rf           = generated_families.rf
}
