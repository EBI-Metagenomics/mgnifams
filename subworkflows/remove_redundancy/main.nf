#!/usr/bin/env nextflow

include { HHSUITE_REFORMAT    } from "${params.moduleDir}/hhsuite/reformat/main.nf"
include { HHSUITE_BUILDHHDB   } from "${params.moduleDir}/hhsuite/buildhhdb/main.nf"
include { POOL_FAMILY_RESULTS } from "${params.moduleDir}/family/main.nf"

workflow REMOVE_REDUNDANCY {
    take:
    fasta_sto_ch
    families_dir

    main:
    a3m_ch = HHSUITE_REFORMAT(fasta_sto_ch, "sto", "a3m").fa
    hh_ch  = HHSUITE_BUILDHHDB(a3m_ch)
    // HHBLITS(...)
    // // non_redundant_families_dir = REMOVE_REDUNDANT(families_dir, hh_ch)
    // pooled_families = POOL_FAMILY_RESULTS(families_dir) // TODO change with non_redundant_families_dir

    // emit:
    // seed_msa_sto = pooled_families.seed_msa_sto
    // msa_sto      = pooled_families.msa_sto
}
