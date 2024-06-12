#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT                     } from "${params.moduleDir}/hhsuite/reformat/main.nf"
include { HHSUITE_BUILDHHDB                    } from "${params.moduleDir}/hhsuite/buildhhdb/main.nf"
include { HHSUITE_HHBLITS                      } from "${params.moduleDir}/hhsuite/hhblits/main.nf"
include { MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID } from "${params.moduleDir}/family/main.nf"
include { POOL_FAMILY_RESULTS                  } from "${params.moduleDir}/family/main.nf"


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
    hh_db_path_ch
    hhr_ch = HHSUITE_HHBLITS(a3m_ch, hh_db_path_ch, db_name).hhr

    mapping = MAP_FIRST_A3M_SEQUENCES_TO_FAMILY_ID(a3m_ch)
    
    // // non_redundant_families_dir = REMOVE_REDUNDANT(families_dir, hh_db_ch)
    // pooled_families = POOL_FAMILY_RESULTS(families_dir) // TODO change with non_redundant_families_dir

    // emit:
    // seed_msa_sto = pooled_families.seed_msa_sto
    // msa_sto      = pooled_families.msa_sto
}
