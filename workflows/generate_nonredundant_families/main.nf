#!/usr/bin/env nextflow

include { GENERATE_FAMILIES_PARALLEL } from "${projectDir}/../../subworkflows/generate_families_parallel/main.nf"
include { MOVE_TO_DIR                } from "${projectDir}/../../modules/family/main.nf"
include { REMOVE_REDUNDANCY          } from "${projectDir}/../../subworkflows/remove_redundancy/main.nf"

workflow {
    families_ch = GENERATE_FAMILIES_PARALLEL( params.clusters_tsv, params.checked_clusters_txt, params.mgnifams_input_fasta_path)

    seed_msa_sto_ch = families_ch.seed_msa_sto.collect()
    msa_sto_ch      = families_ch.msa_sto.collect()
    hmm_ch          = families_ch.hmm.collect()
    rf_ch           = families_ch.rf.collect()
    domtblout_ch    = families_ch.domtblout.collect()
    tsv_ch          = families_ch.tsv.collect()
    discarded_ch    = families_ch.discarded.collect()
    successful_ch   = families_ch.successful.collect()
    converged_ch    = families_ch.converged.collect()
    metadata_ch     = families_ch.metadata.collect()
    logs_ch         = families_ch.logs.collect()

    seed_msa_sto_dir = MOVE_TO_DIR(seed_msa_sto_ch, "seed_msa_sto")
    seed_msa_sto_dir
        .map { filepath ->
            return [ [id:"remove_redundancy"], file(filepath) ]
        }
        .set { seed_msa_sto_dir }
        
    REMOVE_REDUNDANCY(seed_msa_sto_dir, seed_msa_sto_ch, msa_sto_ch, hmm_ch, \
        rf_ch, domtblout_ch, tsv_ch, \
        discarded_ch, successful_ch, converged_ch, \
        metadata_ch, logs_ch)
}
