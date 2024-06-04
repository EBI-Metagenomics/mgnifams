#!/usr/bin/env nextflow

// include { HHSUITE_CREATE_HH_DB } from "${params.moduleDir}/hhsuite/.../main.nf"
// include { REMOVE_REDUNDANT  } from "${params.moduleDir}/family/main.nf"
include { POOL_FAMILY_RESULTS } from "${params.moduleDir}/family/main.nf"

workflow REMOVE_REDUNDANCY {
    take:
    families_dir

    main:
    // hh_models                  = CREATE_HH_DB(families_dir)
    // non_redundant_families_dir = REMOVE_REDUNDANT(hh_models)
    pooled_families = POOL_FAMILY_RESULTS(families_dir) // TODO change with non_redundant_families_dir

    // emit:
    // seed_msa_sto = pooled_families.seed_msa_sto
    // msa_sto      = pooled_families.msa_sto
}
