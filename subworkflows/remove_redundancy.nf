#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT                     } from "../modules/nf-core/hhsuite/reformat/main.nf"
include { HHSUITE_BUILDHHDB                    } from "../modules/nf-core/hhsuite/buildhhdb/main.nf"
include { HHSUITE_HHBLITS                      } from "../modules/nf-core/hhsuite/hhblits/main.nf"
include { COMBINE_HH_RESULTS                   } from "../modules/local/combine_hh_results.nf"
include { MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID } from "../modules/local/map_first_a3m_sequences_to_family_id.nf"
include { POOL_FAM_PROTEINS                    } from "../modules/local/pool_fam_proteins.nf"
include { REMOVE_REDUNDANT_AND_TM              } from "../modules/local/remove_redundant_and_tm.nf"
include { POOL_FAMILY_RESULTS                  } from "../modules/local/family/pool_family_results.nf"

workflow REMOVE_REDUNDANCY {
    take:
    seed_msa_sto_dir
    seed_msa_sto_ch
    msa_sto_ch
    hmm_ch
    rf_ch
    domtblout_ch
    tsv_ch
    discarded_ch
    successful_ch
    converged_ch
    metadata_ch
    logs_ch
    tm_ids_ch
    prob_ids_ch
    rep_fa_ch

    main:
    a3m_ch  = HHSUITE_REFORMAT(seed_msa_sto_dir, "sto", "a3m").fa
    db_name = a3m_ch.map { meta, folderpath ->
        path_str = folderpath.toString()
        parts = path_str.split('/')
        parts[-1]
    }

    hh_db_ch  = HHSUITE_BUILDHHDB(a3m_ch).hh_db
    hh_db_ch
        .map { meta, folderpath ->
            return folderpath
        }
        .set { hh_db_path_ch }
    hhr_ch     = HHSUITE_HHBLITS(a3m_ch, hh_db_path_ch, db_name).hhr
    hhr_all_ch = COMBINE_HH_RESULTS(hhr_ch)

    mapping               = MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID(a3m_ch)
    refined_fam_proteins  = POOL_FAM_PROTEINS(tsv_ch.collect())
    non_redundant         = REMOVE_REDUNDANT_AND_TM(hhr_all_ch, mapping, tm_ids_ch, prob_ids_ch, refined_fam_proteins, rep_fa_ch, params.redundant_threshold, params.similarity_threshold)
    non_redundant_fam_ids = non_redundant.non_redundant_fam_ids
    similarity_edgelist   = non_redundant.similarity_edgelist
    pooled_families       = POOL_FAMILY_RESULTS(seed_msa_sto_ch, \
        msa_sto_ch, hmm_ch, rf_ch, domtblout_ch, tsv_ch, \
        discarded_ch, successful_ch, converged_ch, \
        metadata_ch, logs_ch, non_redundant_fam_ids, similarity_edgelist)

    emit:
    seed_msa_sto = pooled_families.seed_msa_sto
    msa_sto = pooled_families.msa_sto
    hmm = pooled_families.hmm
    rf = pooled_families.rf
    domtblout = pooled_families.domtblout
    tsv = pooled_families.tsv
    discarded = pooled_families.discarded
    successful = pooled_families.successful
    converged = pooled_families.converged
    metadata = pooled_families.metadata
    id_mapping = pooled_families.id_mapping
}
