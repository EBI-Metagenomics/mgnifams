#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT                     } from "${params.moduleDir}/hhsuite/reformat/main.nf"
include { HHSUITE_BUILDHHDB                    } from "${params.moduleDir}/hhsuite/buildhhdb/main.nf"
include { HHSUITE_HHBLITS                      } from "${params.moduleDir}/hhsuite/hhblits/main.nf"
include { COMBINE_HH_RESULTS                   } from "${params.moduleDir}/hhsuite/combine_hh_results.nf"
include { MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID } from "${params.moduleDir}/family/main.nf"
include { REMOVE_REDUNDANT                     } from "${params.moduleDir}/family/main.nf"
include { POOL_FAMILY_RESULTS_FROM_FOLDER      } from "${params.moduleDir}/family/main.nf"

workflow REMOVE_REDUNDANCY {
    take:
    seed_msa_sto_ch
    families_dir

    main:
    a3m_ch  = HHSUITE_REFORMAT(seed_msa_sto_ch, "sto", "a3m").fa
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
    non_redundant         = REMOVE_REDUNDANT(hhr_all_ch, mapping)
    non_redundant_fam_ids = non_redundant.non_redundant_fam_ids
    similarity_edgelist   = non_redundant.similarity_edgelist
    pooled_families       = POOL_FAMILY_RESULTS_FROM_FOLDER(families_dir, non_redundant_fam_ids, similarity_edgelist)
}