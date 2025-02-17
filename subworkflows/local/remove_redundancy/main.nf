#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT                     } from "../../../modules/local/hhsuite/reformat/main"
include { HHSUITE_BUILDHHDB                    } from "../../../modules/local/hhsuite/buildhhdb/main"
include { HHSUITE_HHBLITS                      } from "../../../modules/local/hhsuite/hhblits/main"
include { COMBINE_HH_RESULTS                   } from "../../../modules/local/combine_hh_results/main"
include { MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID } from "../../../modules/local/map_first_a3m_sequences_to_family_id/main"
include { POOL_FAM_PROTEINS                    } from "../../../modules/local/pool_fam_proteins/main"
include { REMOVE_REDUNDANT_AND_TM              } from "../../../modules/local/remove_redundant_and_tm/main"
include { POOL_FAMILY_RESULTS                  } from "../../../modules/local/family/pool_family_results/main"

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
    ch_versions = Channel.empty()

    a3m_ch = HHSUITE_REFORMAT(seed_msa_sto_dir, "sto", "a3m").fa
    ch_versions = ch_versions.mix( HHSUITE_REFORMAT.out.versions )

    db_name = a3m_ch.map { meta, folderpath ->
        def path_str = folderpath.toString()
        def parts = path_str.split('/')
        parts[-1]
    }

    hh_db_ch = HHSUITE_BUILDHHDB(a3m_ch).hh_db
    ch_versions = ch_versions.mix( HHSUITE_BUILDHHDB.out.versions )

    hh_db_ch
        .map { meta, folderpath ->
            return folderpath
        }
        .set { hh_db_path_ch }

    hhr_ch = HHSUITE_HHBLITS(a3m_ch, hh_db_path_ch, db_name).hhr
    ch_versions = ch_versions.mix( HHSUITE_HHBLITS.out.versions )

    hhr_all_ch = COMBINE_HH_RESULTS(hhr_ch)
    // TODO ch_versions = ch_versions.mix( COMBINE_HH_RESULTS.out.versions )

    mapping = MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID(a3m_ch)
    // TODO ch_versions = ch_versions.mix( MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID.out.versions )

    refined_fam_proteins = POOL_FAM_PROTEINS(tsv_ch.collect())
    // TODO ch_versions = ch_versions.mix( POOL_FAM_PROTEINS.out.versions )

    non_redundant = REMOVE_REDUNDANT_AND_TM(hhr_all_ch, mapping, tm_ids_ch, prob_ids_ch, refined_fam_proteins, rep_fa_ch, params.redundant_threshold, params.similarity_threshold)
    // TODO ch_versions = ch_versions.mix( REMOVE_REDUNDANT_AND_TM.out.versions )

    non_redundant_fam_ids = non_redundant.non_redundant_fam_ids
    similarity_edgelist = non_redundant.similarity_edgelist

    pooled_families = POOL_FAMILY_RESULTS(seed_msa_sto_ch, \
        msa_sto_ch, hmm_ch, rf_ch, domtblout_ch, tsv_ch, \
        discarded_ch, successful_ch, converged_ch, \
        metadata_ch, logs_ch, non_redundant_fam_ids, similarity_edgelist, params.starting_id)
    // TODO ch_versions = ch_versions.mix( POOL_FAMILY_RESULTS.out.versions )

    emit:
    versions     = ch_versions
    seed_msa_sto = pooled_families.seed_msa_sto
    msa_sto      = pooled_families.msa_sto
    hmm          = pooled_families.hmm
    rf           = pooled_families.rf
    domtblout    = pooled_families.domtblout
    tsv          = pooled_families.tsv
    discarded    = pooled_families.discarded
    successful   = pooled_families.successful
    converged    = pooled_families.converged
    metadata     = pooled_families.metadata
    id_mapping   = pooled_families.id_mapping
}
